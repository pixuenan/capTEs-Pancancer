#############################################################################
#              1.Identification of noncanonical splicing sites
#############################################################################
## get splicing sites
def get_gtf_intron(gtf, intron_bed):
    # extract intron bed
    cmd = '{gtftk} intronic -i {gtf} --by-transcript -o {out}'.format(
        gtftk=gtftk, gtf=gtf, out=intron_bed)
    os.system(cmd)

## comparing with GENCODE V47
colnames = ['chr','start','end','info','num','strand']
gtfv47_df = pd.read_csv(gtfv47_intron_bed,sep = '\t',header=None,names = colnames)
sample_df = pd.read_csv(sample_intron_bed,sep = '\t',header=None,names = colnames)
gtfv47_df['ss1'] = gtfv47_df['chr']+'_'+gtfv47_df['start'].astype(str)+'_'+gtfv47_df['strand']
gtfv47_df['ss2'] = gtfv47_df['chr']+'_'+gtfv47_df['end'].astype(str)+'_'+gtfv47_df['strand']
sample_df['ss1'] = sample_df['chr']+'_'+sample_df['start'].astype(str)+'_'+sample_df['strand']
sample_df['ss2'] = sample_df['chr']+'_'+sample_df['end'].astype(str)+'_'+sample_df['strand']

## noncanonical splicing sites
gtfv47_SS = list(gtfv47_df['ss1']) + list(gtfv47_df['ss2'])
unique_gtfv47_SS = list(set(gtfv47_SS))
sample_SS = list(sample_df['ss1']) + list(sample_df['ss2'])
unique_sample_SS = list(set(sample_SS))
new_sample_SS = [item for item in unique_sample_SS if item not in set(unique_gtfv47_SS)]

#############################################################################
#        2.Annotation of TE-associated noncanonical splicing sites
#############################################################################
new_sample_SS_df = pd.DataFrame({'Position_site':new_sample_SS})
new_sample_SS_df[['chr','SS_site','strand']] = new_sample_SS_df['Position_site'].str.split('_',expand = True)
new_sample_SS_df['start'] = new_sample_SS_df['SS_site'].astype(int) - 200
new_sample_SS_df['end'] = new_sample_SS_df['SS_site'].astype(int) + 200
new_sample_SS_df_out = new_sample_SS_df[['chr','start','end','strand','SS_site','Position_site']]
new_sample_SS_df_out.to_csv('./Sample_newSS.bed',sep ='\t',index=False,header=False)

## intersect with hg38TE
te_sorted_bed = "./hg38.TE.bed"
SS_bed = './Sample_newSS.bed'

def intersect_te(te_sorted_bed,SS_bed,SS_sorted_bed,out_bed):
    cmd1 = "{bedtools} sort -i {SS_bed} > {SS_sorted_bed}".format(bedtools=bedtools,SS_bed=SS_bed,SS_sorted_bed=SS_sorted_bed)
    cmd2 = "{bedtools} intersect -a {SS_sorted_bed} -b {te_sorted_bed} -wa -wb > {out_bed}".format(bedtools=bedtools,
                                                                                                   SS_sorted_bed=SS_sorted_bed,
                                                                                                   te_sorted_bed=te_sorted_bed,
                                                                                                   out_bed=out_bed)
    os.system(cmd1)
    os.system(cmd2)

## TE-associated noncanonical splicing sites
out_df = pd.read_csv(out_bed,sep='\t',header=None)
out_df_SS = outk_df[[0,1,2,3,4,5]].drop_duplicates()
out_df_SS.to_csv('./Sample_newSS_TE.bed',sep = '\t',header=False,index=False)

#############################################################################
#                      3.Definition of TaPs
#############################################################################
sample_taps_df = sample_df.loc[(sample_df['ss1'].isin(list(out_df_SS[5].unique()))) & (sample_df['ss2'].isin(list(out_df_SS[5].unique())))]









