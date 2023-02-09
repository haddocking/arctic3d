# ARCTIC-3D

[![ci](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Grade)
[![SQAaaS badge shields.io](https://img.shields.io/badge/sqaaas%20software-bronze-e6ae77)](https://api.eu.badgr.io/public/assertions/oAuS52pQTWaC90qMk97hlA "SQAaaS bronze badge achieved")
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8F-orange)](https://fair-software.eu)

[![SQAaaS badge](https://github.com/EOSC-synergy/SQAaaS/raw/master/badges/badges_150x116/badge_software_bronze.png)](https://api.eu.badgr.io/public/assertions/oAuS52pQTWaC90qMk97hlA "SQAaaS bronze badge achieved")

<img src="docs/imgs/arctic3d.png" width="450">

**A**utomatic **R**etrieval and **C**lus**T**ering of **I**nterfaces in Complexes from **3D** structural information

---

## Developing

Check [CONTRIBUTING.md](CONTRIBUTING.md) for more information.

## Installation

### With `conda`

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

And put `blastp` in your `$PATH`.

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
