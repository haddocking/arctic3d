# Benchmarking

This directory can be used to create a benchmark using the [Docking benchmark version 5](https://zlab.umassmed.edu/benchmark/);

```text
python create_benchmark_set Table_BM5.5.xlsx bm5_benchmark.csv
```

This will generate a comma separated table that will serve as the base for the benchmarking. We exclude antibody-antigen complexes and multichain PDBs from the benchmark.

The benchmark file should look like the following:

```csv
complex,receptor,uniprot_receptor,ligand,uniprot_ligand,complex_cat
1AVX_A:B,1QQU_A,P00761,1BA7_B,P01070,EI
1AY7_A:B,1RGH_B,P05798,1A19_B,P11540,EI
1E6E_A:B,1E1N_A,P08165,1CJE_D,P00257,ES
1EWY_A:C,1GJR_A,P21890,1CZP_A,P0A3C8,ES
(...)
```

An [example benchmark file](https://github.com/haddocking/arctic3d/tree/main/benchmark/bm5_benchmark.csv) is already present in the folder.

## run arctic3d on the benchmark

```bash
python3 execute_arctic_bm5.py bm5_benchmark.csv $OUTPUT_DIRECTORY
```

This command runs arctic3d on the full benchmark file `bm5_benchmark.csv`, saving the output in the specified `OUTPUT_DIRECTORY`. More specifically, the target directory will be populated with the results of arctic3d for each unique uniprot ID, together with the arctic3d runs performed in absence of the knowledge about the full complex, that is, excluding the parnter uniprot ID from the search.

The benchmark requires approximately 45 minutes to run on a standard machine.

## HADDOCK and restraints

Output of arctic3d runs can be used as a starting point for some HADDOCK data-driven docking calculations.

In this context, arctic3d output can be employed to generate restraints to drive the docking process, as explained in this [example](https://github.com/haddocking/arctic3d/blob/main/docs/examples.md#arctic3d-restraints).
