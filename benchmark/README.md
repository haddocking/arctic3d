# Benchmarking

## TL:DR

To create a benchmark using the [Docking benchmark version 5](https://zlab.umassmed.edu/benchmark/);

```text
python create_benchmark_set Table_BM5.5.xlsx bm5_benchmark.csv
```

This will generate a comma separated table that will serve as the base for the benchmarking.

```csv
complex,receptor,uniprot_receptor,ligand,uniprot_ligand
1AVX_A:B,1QQU_A,P00761,1BA7_B,P01070
1AY7_A:B,1RGH_B,P05798,1A19_B,P11540
1E6E_A:B,1E1N_A,P08165,1CJE_D,P00257
1EWY_A:C,1GJR_A,P21890,1CZP_A,P0A3C8
(...)
```

## run arctic3d on the benchmark

```bash
python3 execute_arctic_bm5.py bm5_uniprot.csv $OUTPUT_DIRECTORY
```

This command runs arctic3d on the full benchmark file "bm5_uniprot", saving the output in the OUTPUT_DIRECTORY folder.

More specifically, the target directory will be populated with a folder for each complex ID, containing the arctic3d results for the 2 uniprot IDs in play.

The benchmark requires ~ 45 minutes to run on a standard machine.

## Introduction

Coming soon
