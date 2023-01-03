from ast import Str
import streamlit as st
import numpy as np
import pandas as pd
from PIL import Image
import ujson as json
from plotnine import *

# -- Set page config
apptitle = "PhenIEI"

st.set_page_config(page_title=apptitle, page_icon=":genie:", layout="wide")

# -- Set Sidebar
image_pg = Image.open("img/pheniei.png")
st.sidebar.image(image_pg, caption=None, width=100)
st.sidebar.title("PhenIEI")

st.sidebar.header("Exploring knowledge on Inborn Errors of Immunity")

st.sidebar.markdown(
    """
 This webapp browse knowledge on genetic diseases in inborn errors of immunity (IEI).  
 
 If any question or suggestions, please contact: [kevin.yauy@chu-montpellier.fr](kevin.yauy@chu-montpellier.fr)  

 Code source is available in GitHub:
 [https://github.com/kyauy/PhenIEI](https://github.com/kyauy/PhenIEI)

 PhenIEI is an initiative from:
"""
)
image_univ = Image.open("img/logosfacmontpellier.png")
st.sidebar.image(image_univ, caption=None, width=190)

image_chu = Image.open("img/CHU-montpellier.png")
st.sidebar.image(image_chu, caption=None, width=95)


@st.cache(max_entries=50)
def convert_df(df):
    return df.to_csv(sep="\t").encode("utf-8")


@st.cache(allow_output_mutation=True, max_entries=50)
def load_data():
    matrix = pd.read_csv(
        "data/ohe_monarch_iei.tsv.gz",
        sep="\t",
        compression="gzip",
        index_col=0,
    )
    matrix.index = matrix.index.astype(str)
    return matrix


@st.cache(allow_output_mutation=True, max_entries=50)
def load_iei_info():
    matrix = pd.read_csv(
        "data/iei_hpo_2022.tsv",
        sep="\t",
        # compression="gzip",
        index_col=0,
    )
    return matrix


@st.cache(allow_output_mutation=True, max_entries=50)
def add_microdel_to_dict(ncbi_dict):
    new_d = {str(key): str(value) for key, value in ncbi_dict.items()}
    new_d["Del10p13-p14"] = "Del10p13-p14"
    new_d["11q23del"] = "11q23del"
    new_d["14q32"] = "14q32"
    new_d["22q11.2"] = "22q11.2"
    return new_d


@st.cache(allow_output_mutation=True, max_entries=50)
def symbol_to_id_to_dict():
    # from NCBI
    ncbi_df = pd.read_csv("data/Homo_sapiens.gene_info.gz", sep="\t")
    ncbi_df = ncbi_df[ncbi_df["#tax_id"] == 9606]
    ncbi_df_ncbi = ncbi_df.set_index("Symbol")
    ncbi_to_dict_ncbi = ncbi_df_ncbi["GeneID"].to_dict()
    ncbi_df = ncbi_df.set_index("GeneID")
    ncbi_to_dict = ncbi_df["Symbol"].to_dict()
    return add_microdel_to_dict(ncbi_to_dict_ncbi), add_microdel_to_dict(ncbi_to_dict)


@st.cache(
    hash_funcs={"_json.Scanner": hash}, allow_output_mutation=True, max_entries=50
)
def load_hp_ontology():
    with open("data/hpo_obo_202212_iei.json") as json_data:
        data_dict = json.load(json_data)
    return data_dict


@st.cache(allow_output_mutation=True, max_entries=50)
def hpo_description_to_id():
    data_dict = {}
    for key, value in hp_onto.items():
        data_dict[value["name"]] = key
    return data_dict


def get_symbol(gene):
    if gene in symbol.keys():
        return symbol[gene]


def get_hpo_name(hpo):
    names = {}
    if hpo in hp_onto.keys():
        names[hpo] = hp_onto[hpo]["name"]
    return names


def get_hpo_name_only(hpo):
    if hpo in hp_onto.keys():
        return hp_onto[hpo]["name"]
    else:
        return None


def get_hpo_name_list(hpo_list, hp_onto):
    names = {}
    for hpo in hpo_list:
        if hpo in hp_onto.keys():
            names[hpo] = hp_onto[hpo]["name"]
    return dict(sorted(names.items()))


def score(hpo_list, matrix):
    matrix_filter = matrix[hpo_list]
    matrix_filter["sum"] = matrix_filter.sum(axis=1)
    matrix_filter["gene_symbol"] = matrix_filter.index.to_series().apply(get_symbol)
    return matrix_filter.sort_values("sum", ascending=False)


def get_phenotype_specificity(gene_diag, data_patient):
    rank = data_patient.loc[ncbi[gene_diag], "rank"]
    max_rank = data_patient["rank"].max()
    if rank == max_rank:
        return "D - the reported phenotype is NOT consistent with what is expected for the gene/genomic region or not consistent in general."
    elif rank < 23:
        return "A - the reported phenotype is highly specific and relatively unique to the gene (top 23, 50 perc of diagnosis in Simuluk et al. cohort)."
    elif rank < 151:
        return "B - the reported phenotype is consistent with the gene, is highly specific, but not necessarily unique to the gene (top 100, 75 perc of diagnosis in Simuluk et al. cohort)."
    else:
        return "C - the phenotype is reported with limited association with the gene, not highly specific and/or with high genetic heterogeneity."


def get_relatives_list(hpo_list, hp_onto):
    all_list = []
    for hpo in hpo_list:
        all_list.append(hpo)
        if hpo in hp_onto.keys():
            for parent in hp_onto[hpo]["parents"]:
                all_list.append(parent)
            for children in hp_onto[hpo]["childrens"]:
                all_list.append(children)
    return list(set(all_list))


def get_hpo_id(hpo_list):
    hpo_id = []
    for description in hpo_list:
        hpo_id.append(hp_desc_id[description])
    return ",".join(hpo_id)


def add_direct_parents(hpo_list):
    hpo_list_return = []
    for element in hpo_list:
        hpo_list_return.append(element)
        for i in hp_onto[element]["direct_parent"]:
            if hp_onto[i]["distance_to_root"] > 3:
                hpo_list_return.append(i)
    return list(set(hpo_list_return))


hp_onto = load_hp_ontology()
hp_desc_id = hpo_description_to_id()
ncbi, symbol = symbol_to_id_to_dict()
iei_info = load_iei_info()
data = load_data()
gene_list = list(data.index)
symbol_list = [symbol[i] for i in gene_list if i in symbol.keys()]

with st.form("my_form"):
    hpo_raw = st.multiselect(
        "Provide your HPOs",
        list(hp_desc_id.keys()),
        ["Polyarticular arthritis", "Crohn's disease", "Recurrent fever"],
        # [
        #    "Recurrent pneumonia",
        #    "Atopic dermatitis",
        #    "Oligoarthritis",
        #    "Acute hepatitis",
        #    "Recurrent aphthous stomatitis",
        #    "Malar rash",
        #    "Acrocyanosis",
        #    "Increased circulating IgE level",
        # ],
    )

    gene_diag_input = st.multiselect(
        "Optional: provide HGNC gene symbol to be tested",
        options=symbol_list,
        # options=list(ncbi.keys()),
        # default=["STAT3"],
        default=["MEFV"],
        max_selections=1,
    )

    submit_button = st.form_submit_button(
        label="Submit",
    )

if submit_button:
    hpo = get_hpo_id(hpo_raw)
    hpo_list_ini = hpo.strip().split(",")

    if gene_diag_input:
        if gene_diag_input[0] in ncbi.keys():
            gene_diag = gene_diag_input[0]
        else:
            st.write(
                gene_diag_input
                + " gene are not in our database. Please check gene name (need to be in CAPITAL format)."
            )
            gene_diag = None
    else:
        gene_diag = None

    hpo_list_up = []
    for hpo in hpo_list_ini:
        if hpo in ["HP:0000001"]:
            pass
        elif len(hpo) != 10:
            st.write(
                "Incorrect HPO format: "
                + hpo
                + ". Please check (7-digits terms with prefix HP:, and separed by commas)."
            )
            pass
        elif hpo not in data.columns:
            pass
            st.write(hpo + " not available in current database. Please modify.")
        else:
            if data[hpo].astype(bool).sum(axis=0) != 0:
                hpo_list_up.append(hpo)
            else:
                hpo_to_test = hp_onto[hpo]["direct_parent"][0]
                while data[hpo_to_test].astype(bool).sum(
                    axis=0
                ) == 0 and hpo_to_test not in ["HP:0000001"]:
                    hpo_to_test = hp_onto[hpo_to_test]["direct_parent"][0]
                if hpo_to_test in ["HP:0000001"]:
                    st.write(
                        "No gene-HPO associations was found for "
                        + hpo
                        + " and parents."
                    )
                else:
                    hpo_list_up.append(hpo_to_test)
                    st.write(
                        "We replaced: ",
                        hpo,
                        " by ",
                        hp_onto[hpo]["direct_parent"][0],
                        "-",
                        get_hpo_name(hpo_to_test),
                    )
    # hpo_list = list(set(hpo_list_up))
    hpo_list = add_direct_parents(hpo_list_up)

    if hpo_list:
        with st.expander("See HPO inputs"):
            st.write(get_hpo_name_list(hpo_list, hp_onto))

        st.header("Phenotype matching")
        results_sum = score(hpo_list, data).dropna()
        results_sum["matchs"] = results_sum[hpo_list].astype(bool).sum(axis=1)
        results_sum["perc_matchs"] = round(results_sum["matchs"] / len(hpo_list), 2)
        results_sum["score"] = results_sum["matchs"] + results_sum["sum"]
        results_sum["rank"] = (
            results_sum["score"].rank(ascending=False, method="min").astype(int)
        )
        results_sum.loc[results_sum["score"] == 0, "rank"] = 5235
        cols = results_sum.columns.tolist()
        cols = cols[-5:] + sorted(cols[:-5])
        match = results_sum[cols].sort_values(by=["score"], ascending=False)

        st.dataframe(match[match["score"] > 0].drop(columns=["sum"]))
        # st.write(
        #    "Number of genes with at least one match: ", len(match[match["score"] > 0])
        # )
        match_csv = convert_df(match)

        st.download_button(
            "Download matching results",
            match_csv,
            "match.tsv",
            "text/csv",
            key="download-csv-match",
        )

        if gene_diag:
            if ncbi[gene_diag] in results_sum.index:
                st.subheader(
                    "Gene of interest: " + gene_diag,
                )
                p = (
                    ggplot(match, aes("score"))
                    + geom_histogram()
                    + geom_vline(
                        xintercept=results_sum.loc[ncbi[gene_diag], "score"],
                        linetype="dashed",
                        color="red",
                        size=1.5,
                    )
                    + ggtitle("Matching score distribution")
                    + xlab("Gene matching score")
                    + ylab("Number of genes")
                    + theme_bw()
                    + theme(
                        text=element_text(size=12),
                        figure_size=(5, 5),
                        axis_ticks=element_line(colour="black", size=4),
                        axis_line=element_line(colour="black", size=2),
                        axis_text_x=element_text(angle=45, hjust=1),
                        axis_text_y=element_text(angle=60, hjust=1),
                        subplots_adjust={"wspace": 0.1},
                        legend_position=(0.7, 0.35),
                    )
                )
                col1, col2, col3 = st.columns(3)

                with col1:
                    st.pyplot(ggplot.draw(p))

                st.write(
                    "Gene ID rank:",
                    results_sum.loc[ncbi[gene_diag], "rank"],
                    "  |  ",
                    "Gene ID percentage match:",
                    round(results_sum.loc[ncbi[gene_diag], "perc_matchs"] * 100, 2),
                )
                st.write(pd.DataFrame(match.loc[ncbi[gene_diag]]).T)
                st.write(
                    "Gene ID phenotype specificity:",
                    get_phenotype_specificity(gene_diag, results_sum),
                )
            else:
                st.write("Gene ID rank:", " Gene not ranked by PhenoIEI")
            st.subheader(gene_diag + " IUIS descriptions")
            st.write("Disease category")
            st.write(
                iei_info[iei_info["UpdatedGene"] == gene_diag][
                    [
                        "Major category",
                        "Subcategory",
                        "Inheritance",
                    ]
                ]
            )

            st.write("Clinical features")
            st.write(
                iei_info[iei_info["UpdatedGene"] == gene_diag][
                    [
                        "Associated features",
                    ]
                ]
            )
            st.write("Biological features")
            st.write(
                iei_info[iei_info["UpdatedGene"] == gene_diag][
                    [
                        "T cell count",
                        "B cell count",
                        "Immunoglobulin levels",
                        "Neutrophil count",
                        "Other affected cells",
                    ]
                ]
            )
            st.write("Known therapies and pathways")
            st.write(
                iei_info[iei_info["UpdatedGene"] == gene_diag][
                    [
                        "KEGG_drug",
                        "KEGG_pathway",
                    ]
                ]
            )

    else:
        st.write(
            "No HPO terms provided in correct format.",
        )