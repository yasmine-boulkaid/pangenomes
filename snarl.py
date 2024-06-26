import numpy as np
import io
import os
import pandas as pd
from scipy.stats import chi2_contingency
import bdsg
import copy
import re
import matplotlib.pyplot as plt
from tabulate import tabulate
from collections import Counter

my_path = "/home/yboulkaid/Documents/sample_data/pgtest.data/calls/"
all_vcf = os.listdir(my_path)
g0_vcf = []
g1_vcf = []
for i in all_vcf:
    if 'g0' in i:
        g0_vcf.append(i)
    elif 'g1' in i:
        g1_vcf.append(i)

def make_vcf_df(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)), dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, 'FILTER': str, 'INFO': str},
                                    sep ='\t').rename(columns={'#CHROM': 'CHROM'})
def head(df, n = 5, i = 0):
    print(tabulate(df.head(n), headers='keys', tablefmt="mixed_grid"))
def tail(df, n = 5, i = 0):
    print(tabulate(df.tail(n), headers='keys', tablefmt="mixed_grid"))

all_vcf_df = []
for i in all_vcf:
    all_vcf_df.append(make_vcf_df(my_path + i))

g0_vcf_df = []
for i in g0_vcf:
    g0_vcf_df.append(make_vcf_df(my_path + i))

g1_vcf_df = []
for i in g1_vcf:
    g1_vcf_df.append(make_vcf_df(my_path + i))

def chemins(which_vcf_df):
    chemins_possibles = []
    for i in range(len(which_vcf_df["INFO"])):
        text = which_vcf_df["INFO"][i]
        m = re.search('AT=>(.+?);DP', text)
        if m:
            found = m.group(1)
        chemins_possibles.append(found)

    for i in range(len(chemins_possibles)):
        chemins_possibles[i] = chemins_possibles[i].split(',')

    ########################################################################
    chemins_pris = []
    for i in which_vcf_df["SAMPLE"]:
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

    ''' cas simple pour comprendre
    print(chemins_pris[0])
    chemins_pris[0][0] = chemins_possibles[0][chemins_pris[0][0]]
    chemins_pris[0][1] = chemins_possibles[0][chemins_pris[0][1]]
    print(chemins_pris[0])'''

    chemins_possibles = sum(chemins_possibles, [])
    chemins_pris = sum(chemins_pris, [])
    chemins_combines = sum(chemins_combines, [])

    return chemins_possibles, chemins_pris, chemins_combines


def progress_bar(progress, total):
    percent = 100 * (progress / float(total))
    bar = '█' * int(percent) + '-' * (100 - int(percent))
    print(f"\r|{bar}|{percent:.2f}%", end = "\r")

def make_snarl_df(which_vcf_list):
    # chem_fin = chemins_finaux(which_vcf_list)
    df = pd.DataFrame(columns=['snarl index', 'snarl', 'times taken', 'index provisoire'])
    # df['times taken'] = 0 # je sais vraiment pas pourquoi ça marche plus ce machin
    # df = df.assign(times_taken = 0)

    chemins_possibles = [chemins(i)[0] for i in all_vcf_df]
    chemins_possibles = sum(chemins_possibles, [])
    chemins_possibles = list(set(chemins_possibles))
    for i in range(len(chemins_possibles)):
        df.loc[i, 'times taken'] = 0

    df['snarl'] = chemins_possibles


    # Assuming df and which_vcf_list are defined somewhere above this code
    # Also assuming the structure of df and what chemins returns is as expected

    for i in range(len(which_vcf_list)):
        chemins_combines = chemins(which_vcf_list[i])[2]
        combine_count = Counter(chemins_combines)

        # Create a DataFrame from the combine_count dictionary
        count_df = pd.DataFrame(list(combine_count.items()), columns=['snarl', 'count'])

        # Merge this DataFrame with the main DataFrame on the 'snarl' column
        df = df.merge(count_df, on='snarl', how='left', suffixes=('', '_new'))

        # Fill NaN values with 0 in the 'count_new' column
        # df['count_new'] = df['count_new'].fillna(0)

        # Update the 'times taken' column
        df['times taken'] += df['count_new']

        # Drop the 'count_new' column after updating
        df.drop(columns=['count_new'], inplace=True)

    '''for i in range(len(which_vcf_list)):
        chemins_combines = chemins(which_vcf_list[i])[2]
        combine_count = list((x, chemins_combines.count(x)) for x in set(chemins_combines))
        for j in range(len(combine_count)):
            for k in range(len(df['snarl'])):
                if combine_count[j][0] == str(df['snarl'][k]):
                    df.loc[k, "times taken"] += combine_count[j][1]'''
    print('times taken done')

    # fill 'snarl index' column
    for i in range(len(df['snarl index'])):
        if df['snarl'][i][0] == '>':
            S = re.search('>(.+?)>', df['snarl'][i])
            if S:
                s = S.group(1)
            E = re.search('.+>(.*)', df['snarl'][i])
            if E:
                e = E.group(1)
            df.loc[i, "snarl index"] = s + '>' + e
            df.loc[i, "index provisoire"] = int(s)
        else:
            S = re.search('(.+?)>', df['snarl'][i])
            if S:
                s = S.group(1)
            E = re.search('.+>(.*)', df['snarl'][i])
            if E:
                e = E.group(1)
            df.loc[i, "snarl index"] = s + '>' + e
            df.loc[i, "index provisoire"] = int(s)
    print('snarl index done')

    # df.set_index('snarl index', inplace=True, drop=True)
    df.set_index('index provisoire', inplace=True, drop=True)
    df.sort_index(inplace=True)
    df.reset_index()
    df.set_index('snarl index', inplace=True, drop=True)
    print('df done')
    return df


def make_table_contingence(which_snarl_df, which_g0_df, which_g1_df, which_snarl):
    n = which_snarl_df.index.value_counts()[which_snarl]
    chem = []
    for i in range(n):
        chem.append(which_g0_df.loc[which_snarl]['snarl'].iloc[i])
        # thing.append(g0_df.loc[which_snarl]['snarl'][i])
        # corpus_df.loc['it'][1]
    df2 = pd.DataFrame(columns=['g0', 'g1'], index=chem)
    for i in range(len(chem)):
        df2.at[chem[i], 'g1'] = which_g1_df.loc[which_snarl]['times taken'].iloc[i]
        df2.at[chem[i], 'g0'] = which_g0_df.loc[which_snarl]['times taken'].iloc[i]
    return df2


def chi2(table):
    return chi2_contingency(table).pvalue