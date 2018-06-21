#!/usr/bin/python

# ssm_data.txt:
# id: identifier for each SSM. Identifiers must start at s0 and increment, so the first data row will have s0, the second row s1, and so forth.
# gene: any string identifying the variant -- this need not be a gene name. <chr>_<pos> (e.g., 2_234577) works well.
# a: number of reference-allele reads at the variant locus.
# d: total number of reads at the variant locus.
# mu_r: fraction of expected reference allele sampling from the reference population. E.g., if the tumor has an A->T somatic mutation at the locus, the genotype of the reference population should be AA. Thus, mu_r should be 1 - (sequencing error rate). Given the 0.001 error rate in Illumina sequencing, setting this column to 0.999 works well.
# mu_v: fraction of expected reference allele sampling from variant population. Suppose an A->T somatic mutation occurred at the locus. mu_v always uses normal ploidy (i.e., the copy number in non-CNV regions). As humans are diploid, copy number will thus always be 2. So, the variant population genotype should be AT, meaning we will observe the reference allele with frequency 0.5 - (sequencing error rate). Given the 0.001 error rate in Illumina sequencing, setting this column to 0.499 works well.

import sys
import pandas as pd

if __name__ == "__main__":

    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python %s [data_file]\n" % sys.argv[0])
        sys.exit(1)

    filename = sys.argv[1]
    COVERAGE=1000000

    with open(filename) as f:
        f.readline()
        k = int(f.readline().split()[0])
        n = int(f.readline().split()[0])

    df = pd.read_table(filename, skiprows=3)
    sys.stdout.write("\t".join(['id', 'gene', 'a', 'd', 'mu_r', 'mu_v']) + "\n")
    for c in range(n):
        identifier = "s%d" % c
        gene = df[df['character_index'] == c].iloc[0]['character_label']

        sys.stdout.write("\t".join(map(str, [identifier, gene])) + "\t")

        for p in range(k):
            f = df[(df['character_index'] == c) & (df['#sample_index'] == p)]['f-']
            a = COVERAGE - int(f / 2 * COVERAGE)
            if p > 0:
                sys.stdout.write(",")
            sys.stdout.write(str(a))

        sys.stdout.write("\t" + ",".join([str(COVERAGE) for p in range(k)]))
        sys.stdout.write("\t1\t0.5\n")
