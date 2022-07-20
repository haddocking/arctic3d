# ARCTIC-3D

[![python lint](https://github.com/haddocking/arctic3d/actions/workflows/.lint.yml/badge.svg)](https://github.com/haddocking/arctic3d/actions/workflows/.lint.yml)
[![unittests](https://github.com/haddocking/arctic3d/actions/workflows/unittests.yml/badge.svg)](https://github.com/haddocking/arctic3d/actions/workflows/unittests.yml)

<img src="docs/imgs/arctic3d.png" width="450">

**A**utomatic **R**etrieval and **C**lus**T**ering of **I**nterfaces in Complexes from **3D** structural information

---

## Developing

Check [DEVELOPING.md](DEVELOPING.md) for more information.

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

## Example

```bash
# the default input is an uniprotID
arctic3d P00760
# or to identify the uniprotID number remotelly
arctic3d example/1ppe_E.fasta
# or to identify the uniprotID number locally
arctic3d example/1ppe_E.fasta --db db/swissprot
```

---
