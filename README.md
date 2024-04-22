# ARCTIC-3D

![PyPI - License](https://img.shields.io/pypi/l/arctic3d)
![PyPI - Status](https://img.shields.io/pypi/status/arctic3d)
![PyPI - Version](https://img.shields.io/pypi/v/arctic3d)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/arctic3d)

[![ci](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Grade)

[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-bronze-e6ae77)](https://api.eu.badgr.io/public/assertions/oAuS52pQTWaC90qMk97hlA "SQAaaS bronze badge achieved")
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)


[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_bronze.png)](https://api.eu.badgr.io/public/assertions/oAuS52pQTWaC90qMk97hlA "SQAaaS bronze badge achieved")

<img src="https://raw.githubusercontent.com/haddocking/arctic3d/main/docs/imgs/arctic3d.png" width="450">

**A**utomatic **R**etrieval and **C**lus**T**ering of **I**nterfaces in Complexes from **3D** structural information

## WEB SERVER

ARCTIC-3D is available at this webserver https://wenmr.science.uu.nl/arctic3d/

## ARCTIC-3D: all you want to know about protein-specific interfaces

ARCTIC-3D is a software for data-mining and clustering of protein interface information. It allows you to retrieve all the existing interface information for your desired protein from the PDBE graph database (https://www.ebi.ac.uk/pdbe/pdbe-kb/), grouping similar interfaces in interacting surfaces.

The software first checks your input (a uniprot ID, a FASTA file, or a PDB file), and then retrieves the existing interaction data from the [graph API](https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/). Such interfaces are projected on a selected PDB structure and their dissimilarity is calculated, thus allowing for the application of a hierarchical clustering algorithm.

In output you will see how your favourite protein can display different binding surfaces, each one characterised by few residues that are always present (_hotspots_) and other amino acids which are at the interface only from time to time.

## Developing

Check [CONTRIBUTING.md](CONTRIBUTING.md) for more information.

## Installation

### With `conda`

Clone the repository on your computer and navigate to it

```bash
git clone git@github.com:haddocking/arctic3d.git
cd arctic3d
```

Here you can create the arctic3d environment:

```bash
conda create -n arctic3d python=3.10
conda activate arctic3d
pip install .
arctic3d -h
```

## To run BLAST locally

```bash
bash install_blast_deps.sh
```

And put `blastp` in your `$PATH` by adding the following line to your `.bashrc` or `.bash_profile` file:

```bash
export PATH="PATH_TO_YOUR_ARCTIC3D_INSTALLATION/src/ncbi-blast-2.15.0+/bin:$PATH"
```

## Example usage

Please refer to the [examples](docs/examples.md) documentation page.

## Detailed documentation

In order to generate a detailed html documentation please execute these commands

```text
pip install myst_parser
pip install chardet
conda install sphinx
sphinx-build -E docs ./arctic3d-docs
```

Then you can open the file `arctic3d-docs/index.html`, which contains all the necessary documentation.

## Citing us

If you used ARCTIC-3D in your work please cite the following publication:

**Marco Giulini, Rodrigo V. Honorato, Jesús L. Rivera, and Alexandre MJJ Bonvin**: "ARCTIC-3D: automatic retrieval and clustering of interfaces in complexes from 3D structural information." Communications Biology 7, no. 1 (2024): 49. (www.nature.com/articles/s42003-023-05718-w)
