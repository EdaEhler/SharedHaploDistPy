# SharedHaploDistPy
Compute Shared Haplogroup Distance (SHD) and associated p-values.
This distance method to be used in population genetics with frequencies of mitochondrial haplogroups were published by Neparáczki et al. (2018) (https://doi.org/10.1371/journal.pone.0205920)

### SHD computation examples:
this will create 2 new output files, (i) one with a matrix of computed SHD  pair-wise distances (`SHDoutput_dist.csv`), (ii) the second with a matrix of their p-values created by 10000 permutation steps (`SHDoutput_pvals.csv`):

```
python SharedHaploDistPy.py -i SHD_example_input.csv -s ';' -p 10000 -o SHDoutput
```
