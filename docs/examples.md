# Execution Examples

ARCTIC3D can be executed in several ways:

1. [standard-uniprot-input](#standard-uniprot-input)
1. [standard-fasta-input](#standard-fasta-input)
1. [standard-pdb-input](#standard-pdb-input)
1. [uniprot-full-input](#uniprot-full-input)
1. [uniprot-ligand-input](#uniprot-ligand-input)
1. [arctic3d-resclust](#arctic3d-resclust)
1. [arctic3d-localise](#arctic3d-localise)
1. [arctic3d-restraints](#arctic3d-restraints)

## standard-uniprot-input

In this example a uniprot ID is provided as input to arctic3d.

```bash
arctic3d P00760
```

This is the basic ARCTIC3D scenario: here the program performs the following tasks sequentially:

1. retrieval of all the partners interacting with the selected uniprot ID;
2. gathering of the corresponding interfaces
3. structural projection of the interfaces on the [best available pdb](https://www.ebi.ac.uk/pdbe/api/doc/sifts.html)
4. interface clustering using the irregular spin glass theory

Step 3 can be skipped by specifying the pdb that must be retrieved:

```bash
arctic3d P00760 --pdb_to_use=4xoj
```

This typically allows to considerably reduce the execution time of the program.

One or more uniprot IDs can be excluded from set of interacting partners. As an example, this may be important if you want to exclude homodimer interfaces from the search:

```bash
arctic3d P00760 --out_uniprot=P00760,P00974
```

Following the same logic, one or more pdb IDs can be excluded from the search:

```bash
arctic3d P00760 --out_pdb=4xoj,6sy3
```

## standard-fasta-input

Here a single sequence (in FASTA format) is given to ARCTIC3D. The program blasts the sequence against the local version of blast to get the corresponding uniprot ID.

```bash
arctic3d example/1ppe_E.fasta --db db/swissprot
```

Once the uniprot ID is retrieved the workflow proceeds as explained in [standard-uniprot-input](standard-uniprot-input).

To run blast remotely simply remove the db reference:

```bash
arctic3d example/1ppe_E.fasta
```

## standard-pdb-input

It is also possible to provide ARCTIC3D with a set of pre-calculated interfaces. These might be extracted with any computational or experimental methodology. In this case an input pdb is mandatory:

```bash
arctic3d example/1ppe_E.pdb --interface_file example/1ppe_E_example_interfaces.txt
```

One may want to run ARCTIC3D on a custom pdb, without providing an interface file. Here the program extracts the pdb sequence and finds the corresponding uniprot ID. The input file is then renumbered according to the canonical numbering. Finally, the workflow showcased in [standard-uniprot-input](standard-uniprot-input) can be executed.

## uniprot-full-input

In the PDB it might happen that the same protein partners interact with a different set of residues, depending on a variety of factors. As an example, proteins forming homotrimers form at least two different interfaces. By default, ARCTIC3D does not consider this multimeric information, but the `full` option can be utilized to separate protein-protein interfaces at the pdb-chain level.

```bash
arctic3d P00760 --full=True
```

## uniprot-ligand-input

Small-molecule interaction information is not retrieved by default by ARCTIC3D. By setting the `ligand` option to `yes` the retrieval can be limited to ligands:

```bash
arctic3d P00760 --ligand=yes
```

Such information can be combined with standard interface information by setting the `ligand` option to `both`:

```bash
arctic3d P00760 --ligand=both
```

**Disclaimer**: typically, different small molecules of the same type (ions in particular) cannot be distinguished according to the [graph-api](https://www.ebi.ac.uk/pdbe/pdbe-kb/api). Therefore results have to be interpreted carefully.

## arctic3d-resclust

It is also possible to use ARCTIC3D to cluster separate residues on a protein structure:

```bash
arctic3d-resclust example/1ppe_E.pdb --residue_list 49,50,51,100,101,102 --threshold=12.0 --chain=E --linkage=average
```

Here each residue is treated as an independent entity and the standard CA-CA distance matrix between the selected amino acids is clustered according to the `threshold` value and the `linkage` criterion. This can be useful if one wants to automatically separate groups of residues on a pdb structure.

The user can change the `linkage` strategy employed to create the dendrogram by choosing between one of the keywords specified [on the scipy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html).

## arctic3d-localise

Did you run ARCTIC3D on your favourite protein and wonder how the different interface clusters depend on the subcellular localisation of the partners? Try running `arctic3d-localise` on ARCTIC3D output files:

```bash
arctic3d-localise ${arctic3d_rundir}/clustered_interfaces.out
```

Here `arctic3d_rundir` is the name of the output folder that you want to analyse.

## arctic3d-restraints

It is possible to use the results of two ARCTIC3D runs to generate HADDOCK-specific restraints to inform a data-driven docking process:

```bash
arctic3d-restraints --r1 ${first_arctic3d_rundir} --r2 ${second_arctic3d_rundir} --prob_threshold=0.5 --ch1=A --ch2=B
```

This works by extracting the residues that are more likely to occur in the clusters identified by ARCTIC3D, according to the selected probability threshold (0.3 by default). In fact, it does not always make sense to include all the residues found in a cluster. Parameters `ch1` and `ch2` define the chain ID that will be used in the restraint files for each interacting partner.

The output of this command is a new folder (`arctic3-restraints` by default) that contains all the generated restraint .tbl files, one for each pair of interface clusters. For example, if there are three clusters in `first_arctic3d_rundir` and four in `second_arctic3d_rundir`, you will obtain 12 different restraint files. The .tbl files are also zipped in a tgz archive (`ambig.tbl.tgz`), which can be directly [used within HADDOCK3](https://github.com/haddocking/haddock3/blob/main/docs/examples.md#docking-multiple-ambig).
