
from variables import *

vcf1 = make_vcf_df(my_path + all_vcf[0])

all_df = make_snarl_df(all_vcf_df)
g0_df  = make_snarl_df(g0_vcf_df)
g1_df  = make_snarl_df(g1_vcf_df)
head(all_df)
head(g0_df)
head(g1_df)



