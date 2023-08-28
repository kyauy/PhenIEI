---
title: PhenIEI
emoji: fire
sdk: streamlit
sdk_version: 1.21.0 
app_file: pheniei_app.py
pinned: true
---

# PhenIEI
Exploring knowledge on Inborn Errors of Immunity

![](img/pheniei.png)


Contact : [kevin.yauy@chu-montpellier.fr](mailto:kevin.yauy@chu-montpellier.fr)


## Introduction

Through next-generation sequencing adoption, the clinical and genetic spectrum of Inborn Errors of Immunity (IEI) continuously expands with overlapping phenotypes across its ten major categories of diseases. Increasing phenotypic and genetic variability of known IEI is a critical challenge for IEI diagnosis. 

We extracted knowledge on genetic diseases in inborn errors of immunity (IEI) based on [2022 IUIS classification](https://link.springer.com/article/10.1007/s10875-022-01289-3), [GenIA database](https://www.geniadb.net/) and [HPO](https://hpo.jax.org/app/). 

Based on the gathered knowledge, we developed a clinical decision support system called *PhenIEI* to match patients’ symptoms to IEI knowledge, where genes are ranked according to the number of matched symptoms.


## Run the tool

A webapp is accessible at [https://pheniei.streamlit.app/](https://pheniei.streamlit.app/), **please try it !**

It's a streamlit application, where code is accessible in ̀`pheniei_app.py` file. 

To install on your local machine, you need `poetry` package manager and launch in the folder:
```
poetry install
```

To make it run in your local computer:
```
poetry shell
streamlit run pheniei_app.py
```

Using requirement ?
```
poetry export --without-hashes --format=requirements.txt > requirements.txt
```

--------------------------------------------------------------------------------
 
*PhenIEI is an initiative from:*
- [Université de Montpellier - LIRMM](https://www.igmm.cnrs.fr/service/equipe-igmm-lirmm-imag-regulations-genomiques-computationnelles/)
- [CHU de Montpellier - Unité Médicale des Maladies Rares et Auto-Inflammatoires](https://umai-montpellier.fr/)

*PhenIEI is in partnership with:*
- [GenIA database](https://www.geniadb.net/)