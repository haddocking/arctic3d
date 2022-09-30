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
# one or more uniprot IDs can be excluded from the interface calls
arctic3d P00760 --out_uniprot=P00760,P00974
# one or more pdb IDs can be excluded from the interface calls
arctic3d P00760 --out_pdb=4xoj,6sy3
# the pdb that must be retrieved can be specified
arctic3d P00760 --pdb_to_use=4xoj
# the uniprotID can be identifed from the sequence, either remotely
arctic3d example/1ppe_E.fasta
# or locally
arctic3d example/1ppe_E.fasta --db db/swissprot
# it is also possible to provide a pdb file with a set of interfaces
arctic3d example/1ppe_E.pdb --interface_file example/1ppe_E_example_interfaces.txt
```

---

## Residue-based clustering

It is also possible to clusters separate residues with the following command:
```bash
arctic3d example/1ppe_E.pdb ./example/1ppe_E.pdb --residue_list 49,50,51,100,101,102 --threshold=12.0 --segid=E
```
Here each residue is treated as an independent entity and the standard CA-CA distance matrix between the selected amino acids is clustered according to the `threshold` value.
This can be useful if one wants to automatically separate groups of residues on a pdb structure.
