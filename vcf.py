import pandas as pd
import io
import os
import re
import copy
from tabulate import tabulate


class Vcf:
    path = "/home/yboulkaid/Documents/sample_data/pgtest.data/calls/"
    def __init__(self):
        self.df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

    def make_table(self, path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        self.df = pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
                                        sep ='\t').rename(columns={'#CHROM': 'CHROM'})

    def head(self, n = 5, i = 0):
        print(tabulate(self.df.head(n), headers='keys', tablefmt="mixed_grid"))

    def tail(self, n = 5, i = 0):
        print(tabulate(self.df.tail(n), headers='keys', tablefmt="mixed_grid"))

    def chemins(self):
        chemins_possibles = []
        for i in range(len(self.df["INFO"])):
            text = self.df["INFO"][i]
            m = re.search('AT=>(.+?);DP', text)
            if m:
                found = m.group(1)
            chemins_possibles.append(found)

        for i in range(len(chemins_possibles)):
            chemins_possibles[i] = chemins_possibles[i].split(',')

        ########################################################################
        chemins_pris = []
        for i in self.df["SAMPLE"]:
            found = i[0:3]
            chemins_pris.append(found)

        for i in range(len(chemins_pris)):
            chemins_pris[i] = chemins_pris[i].split('/')
            for j in range(len(chemins_pris[i])):
                chemins_pris[i][j] = int(chemins_pris[i][j])

        ########################################################################
        chemins_combines = copy.deepcopy(chemins_pris)

        for i in range(len(chemins_pris)):
            chemins_combines[i][0] = chemins_possibles[i][chemins_pris[i][0]]
            chemins_combines[i][1] = chemins_possibles[i][chemins_pris[i][1]]

        chemins_possibles = sum(chemins_possibles, [])
        chemins_pris = sum(chemins_pris, [])
        chemins_combines = sum(chemins_combines, [])

        return chemins_possibles, chemins_pris, chemins_combines


all_vcf = os.listdir(Vcf.path)
g0_vcf = []
g1_vcf = []
for i in all_vcf:
    if 'g0' in i:
        g0_vcf.append(i)
    elif 'g1' in i:
        g1_vcf.append(i)

all_vcf_df = []
for i in all_vcf:
    vcf = Vcf()
    vcf.make_table(Vcf.path + i)
    all_vcf_df.append(vcf.df)

g0_vcf_df = []
for i in g0_vcf:
    vcf = Vcf()
    vcf.make_table(Vcf.path + i)
    g0_vcf_df.append(vcf.df)

g1_vcf_df = []
for i in g1_vcf:
    vcf = Vcf()
    vcf.make_table(Vcf.path + i)
    g1_vcf_df.append(vcf.df)


