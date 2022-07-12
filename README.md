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

## Example

```bash
arctic3d example/1ppe_E.fasta
```

---
