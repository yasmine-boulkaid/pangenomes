from scipy.stats import chi2_contingency

class Table:
    def __init__(self, snarl):
        self.snarl = snarl

    def make_table(self):
        n = self.snarl

    def table_contingence(which_snarl):
        n = snarl_df.index.value_counts()[which_snarl]
        chem = []
        for i in range(n):
            chem.append(g0_df.loc[which_snarl]['snarl'].iloc[i])
            # thing.append(g0_df.loc[which_snarl]['snarl'][i])
            # corpus_df.loc['it'][1]
        df2 = pd.DataFrame(columns=['g0', 'g1'], index=chem)
        for i in range(len(chem)):
            df2.at[chem[i], 'g1'] = g1_df.loc[which_snarl]['times taken'].iloc[i]
            df2.at[chem[i], 'g0'] = g0_df.loc[which_snarl]['times taken'].iloc[i]
        return df2

    def chi2(self):
        return chi2_contingency(self).pvalue