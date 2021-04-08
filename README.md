---
output:
  pdf_document: default
  html_document: default
---

# ctnet
An `R` package to aid researchers analyse experience sampling data using a Continuous-Time Dynamical Network approach. The package takes output from CT-VAR models fit using the `ctsem` package and processes the output, allowing researchers to obtain estimated path-specific effects, centrality measures, the uncertainty around those measures, and to simulate press and pulse interventions based on the estimated model.  

## Background
This repository contains an `R` package used by Ryan \& Hamaker (in press)[[PsyArXiv]](https://psyarxiv.com/ryg69/) to aid researchers in conducting dynamical network analysis using continuous-time statistical models.
(here referred to as weighted DAG structures). 


## Installation
The current version of this package can be installed directly from github using
```r
devtools::install_github("ryanoisin/ctnet")
```

## Usage
The main functions take output from `ctsem` models as input. 

## Contact Details

For more details please contact **o.ryan@uu.nl**
