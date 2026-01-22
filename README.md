# ARCTIC-3D

![PyPI - License](https://img.shields.io/pypi/l/arctic3d)
![PyPI - Status](https://img.shields.io/pypi/status/arctic3d)
![PyPI - Version](https://img.shields.io/pypi/v/arctic3d)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/arctic3d)

[![ci](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml/badge.svg)](https://github.com/haddocking/arctic3d/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/dc788367452c47928e30f2f1f481d7e4)](https://www.codacy.com/gh/haddocking/arctic3d/dashboard?utm_source=github.com&utm_medium=referral&utm_content=haddocking/arctic3d&utm_campaign=Badge_Grade)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10839520.svg)](https://doi.org/10.5281/zenodo.10839520)
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

## Result interpretation

After running ARCTIC-3D, results are stored in the output directory (default: `arctic3d-{uniprot_id}/`). Below is an explanation of each output file and how to interpret them.

### Output files

| File | Description |
|------|-------------|
| `arctic3d.log` | Log file with execution details and warnings |
| `input_data/` | Directory containing copies of input files |
| `{pdb_id}_updated.cif` | Structure file downloaded from PDBe (mmCIF format) |
| `{pdb_id}-{chain}.pdb` | Cleaned PDB structure used for analysis (renumbered to UniProt numbering) |
| `retrieved_interfaces.out` | All interfaces retrieved from PDBe, listing partner IDs and their residue lists |
| `interface_matrix.txt` | Pairwise dissimilarity values between all interfaces (used for clustering) |
| `dendrogram_{linkage}.png` | Hierarchical clustering dendrogram (e.g., `dendrogram_average.png`) |
| `clustered_interfaces.out` | Interfaces grouped into clusters (binding surfaces) |
| `clustered_residues.out` | Residues belonging to each cluster |
| `clustered_residues_probs.out` | Residues ranked by probability within each cluster |
| `{pdb_id}-{chain}_cl{N}.pdb` | PDB structure for cluster N with probabilities encoded in B-factor column |
| `sequence_probability.html` | Interactive bar plot of per-residue probabilities |
| `sequence_probability.json` | JSON data for the interactive plot |

### Understanding the clustering

ARCTIC-3D groups similar interfaces into **binding surfaces** (clusters). Two interfaces are considered similar when they overlap spatially on the protein surface. The dissimilarity is measured using the squared sine of the angle between interface vectors in a Hilbert space representation - values close to 0 indicate overlapping interfaces, while values close to 1 indicate completely distinct regions.

The `interface_matrix.txt` file contains the pairwise dissimilarity values in the format:
```
interface1 interface2 dissimilarity_value
```

### Interpreting residue probabilities

The **probability** (or "contact probability score") represents **the fraction of interfaces within a cluster where a residue is observed**. It is calculated independently for each cluster:

```
probability = (number of interfaces containing the residue) / (total interfaces in cluster)
```

For each cluster, residues are assigned a probability value between 0 and 1:

- **Probability = 1.0**: The residue appears in every interface within the cluster (a "hotspot" residue)
- **Probability = 0.5**: The residue appears in half of the cluster's interfaces
- **Probability close to 0**: The residue rarely appears at this binding surface

**Important**: Probabilities do NOT sum to 1.0 across clusters for a given residue. A residue can have high probability in multiple clusters if it participates in different binding surfaces. For example, a residue with probability 0.8 in cluster 1 and 0.6 in cluster 2 means it appears in 80% of cluster 1's interfaces and 60% of cluster 2's interfaces.

The `clustered_residues_probs.out` file lists residues ranked by probability:
```
Cluster 1 : 15 residues
rank    resid   resname probability
1       42      ALA     1.000
2       45      GLU     0.875
...
```

### Visualizing probabilities in PDB files

The output PDB files (`{pdb_id}-{chain}_cl{N}.pdb`) encode probabilities in the B-factor column:

- Cluster residues: `B = 50 × (1 + probability)`, ranging from 50 (probability=0) to 100 (probability=1)
- Non-cluster residues: `B = 0`

This allows visualization in molecular viewers (PyMOL, ChimeraX, etc.) using a color spectrum where high B-factors (red) indicate hotspot residues and low values (blue) indicate residues not involved in that binding surface.

### Interpreting the dendrogram

The dendrogram (`dendrogram_average.png`) shows the hierarchical relationship between all retrieved interfaces. The x-axis represents the dissimilarity between interfaces or groups. Interfaces that merge below the threshold (default: 0.866, corresponding to a 60° angle) form a single binding surface. The threshold can be adjusted with `--threshold` to obtain finer or coarser clustering.

### Using the interactive plot

Open `sequence_probability.html` in a web browser to explore per-residue binding probabilities. Each cluster is shown as a separate colored bar series, allowing you to identify which residues are involved in which binding surface and compare hotspots across different clusters.

## Citing us

If you used ARCTIC-3D in your work please cite the following publication:

**Marco Giulini, Rodrigo V. Honorato, Jesús L. Rivera, and Alexandre MJJ Bonvin**: "ARCTIC-3D: automatic retrieval and clustering of interfaces in complexes from 3D structural information." Communications Biology 7, no. 1 (2024): 49. (www.nature.com/articles/s42003-023-05718-w)
