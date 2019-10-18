# command line interface (using argparse) use for SHD computing

# IMPORTS
import argparse
import pandas as pd
import numpy as np

from random import shuffle # for pvalues computation/permutation
from itertools import combinations # cycling through population pairs
from scipy.spatial.distance import pdist, squareform # compute and parse SHD matrix


# FUNCTIONS
def compute_shd(pop1, pop2):
    """
    Computes the Shared Haplogroup Distance between 2 numpy/pandas vectors.
    Haplotype frequencies needed.
    """
    assert len(pop1) == len(pop2), "Length of pop1 and pop2 vector do not match!"
    
    # 1/2 * sum( |xi - yi| )
    shd = sum(abs(pop1 - pop2))/2
    
    return shd




def permute_shd_pval(pop1, pop2, n_permutations=1000, counts=None, dist=None):
    """
    Computes the statistical significance of Shared Haplotype Distance (SHD). Same approach as in
    Arlequin 3.5 to obtaining Fst p-values.
    The null distribution of pairwise SHD values under the null hypothesis (no difference between
    the two populations) is obtained by permuting haplotypes between populations. The p-value of
    the test is the proportion of permutations leading to a SHD value larger or equal to the
    observed one.
    pop1, pop2 - names of populations, used as index in counts pandas.DataFrame
    n_permutations - how many permutations to create the null distribution
    counts - pd.DataFrame, first column is n, the rest are mt hg counts
    dist - the original SHD value to be tested
    """

    count = 0  # number of distances higher than 'dist' after permutation

    pop1_n = counts.loc[pop1, "n"]
    pop2_n = counts.loc[pop2, "n"]

    # order of columns for distance permuatations computation
    hgs = ["A", "B", "D", "F", "G", "M", "C", "Z", "N*", "N1a",
           "I", "W", "X", "R", "HV", "V", "H", "T", "J", "U",
           "U2", "U3", "U4", "U5a", "U5b", "U8", "K"]

    ser1 = counts.loc[pop1][1:]
    ser2 = counts.loc[pop2][1:]

    # for each population make a list with so many hg symbols (i.e. 'A')
    # that each population has counts -> double list comprehension
    pool = [y for x in zip(ser1.index, ser1) for y in [x[0]] * x[1]] + [y for x in zip(ser2.index, ser2) for y in [x[0]] * x[1]]

    # main permutation loop
    for n in range(n_permutations):
        shuffle(pool)  # mix the two population lists
        new_p1 = pool[:pop1_n]  # split according to population sizes (original pop1)
        new_p2 = pool[pop1_n:]  # split according to population sizes (original pop2)
        # print(new_p1, new_p2)

        c1 = [new_p1.count(i) for i in hgs]  # count hg symbols in new populations
        c2 = [new_p2.count(i) for i in hgs]

        # compute_shd with the new counts (turn them to frequencies first (np.array))
        new_dist = compute_shd(np.array([i/pop1_n for i in c1]), np.array([i/pop2_n for i in c2]))
        count += (new_dist >= dist)  # add 1 to count if True

    return count / n_permutations  # p-value


def command_line_interface():
    # command line interface

    intro = """\
Compute Shared Haplogroup Distance (SHD) and associated p-values
================================================================
    """


    examples = """\
SHD computation examples:
-------------------------

# this will create 2 new files with computed distance
# and their p-values created by 10000 permutation steps:
# SHDoutput_pvals.csv and SHDoutput_dist.csv

python SharedHaploDistPy.py -i SHD_example_input.csv -s ';' -p 1000 -o SHDoutput

    """


    # initiate parser, more arg groups created
    parser = argparse.ArgumentParser(description=intro,
                                     epilog=examples,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input",
                        help="input file, delimited table, first columns are population, n, \
                        followed by mitochondrial haplogroups (HG) counts \
                        (how many samples from the given population \
                        are from each  HG.", required=True)

    parser.add_argument("-o", "--output", help="output file prefix, if not present, \
                        input file name will be parsed for ouput file name",
                        default=None)
    parser.add_argument("-s", "--separator", help="separator character for the \
                        input file (i.e. ';', ',', '\\t')", default=";")
    parser.add_argument("-p", "--permutations", help="number of permutation steps \
                        used to compute the p-values.", default=1000, type=int)

    # parsing of the sys.argv
    args = parser.parse_args()

    return args.input, args.output, args.separator, args.permutations

    # draw_heatmap(args.input, pic_name=args.output, title=args.title,
    #             y_label=args.ylabel, cutoff=args.cutoff,
    #             marker_x=args.markerx, marker_label=args.markerlabel,
    #             dpi=args.dpi)  # add more arguments

########################################
# compute SHD matrix + p-values matrix #
########################################

# data-structures needed:
# -----------------------
# counts - matrix of hg counts (user input)
# data - matrix of hg freqs without the 'n' column (parse from counts)

infile, outfile, separator, permutations = command_line_interface()

if not outfile:
    # someInput.Name.csv --> "someInput.Name"
    outfile = infile.rsplit(".", 1)[0]

# read input, first column will be an index
counts = pd.read_csv(infile, sep=separator, index_col=0)

print("\n\n")
print("Computing Shared Haplogroup Distance (SHD) and associated p-values")
print("==================================================================")
print()
print("Read input:")
print("-----------")
print(counts)
print("\n\n")

# do not use 'n' column
data = counts.apply(lambda x: x / x["n"], axis=1)
data = data.iloc[:, 1:]

# print(data)

shd = pdist(data.values, metric=compute_shd)  # vectorized, very fast ~1 sec
shd_matrix = squareform(shd)
shd_distances = pd.DataFrame(shd_matrix,
                             index=counts.index,
                             columns=counts.index.values)
pvals_df = pd.DataFrame(index=counts.index, columns=counts.index.values)

print("SHD distances:")
print("--------------")
print(shd_distances)
print("\n\n")

counter = 0

for ind, col in combinations(counts.index, 2):
    pvals_df.loc[ind, col] = permute_shd_pval(ind, col,
                                              n_permutations=permutations,
                                              dist=shd_distances.loc[ind, col],
                                              counts=counts)
    counter += 1

print("SHD p-values (upper triangle):")
print("------------------------------")
print(pvals_df)
print("\n\n")

# output the results
pvals_df.to_csv(outfile + "_pvals.csv", sep=separator)
shd_distances.to_csv(outfile + "_dist.csv", sep=separator)
