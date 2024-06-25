import os
from vcf import Vcf, all_vcf_df, g0_vcf_df, g1_vcf_df
import table

vcf1 = Vcf()
vcf1.make_table("/home/yboulkaid/Documents/sample_data/pgtest.data/calls/samp_g0_0.vcf")

chemins_possibles1, chemins_pris1, chemins_combines1 = vcf1.chemins()
print('========== possibles ==========')
print(chemins_possibles1[0:5])
print('============ pris =============')
print(chemins_pris1[0:5])
print('========== combinés ===========')
print(chemins_combines1[0:5])
