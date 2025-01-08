import pandas as pd
from numpy import *
import pickle, os, sys
#with open('Human_hg38_10kb_windows.pkl','rb') as a:
#    W = pickle.load(a)

keyword = sys.argv[1]
group = sys.argv[2]
outfile = 'merged_table_no_correction_'+group+'_10kb.npz'

with open('excluded_bands.pkl','rb') as a:
    exclude_bands = pickle.load(a)

with open('Human_hg38_auto_10kb_windows.pkl','rb') as b:
    W = pickle.load(b)

exclude_keys = []
for r in W:
    if W[r] in exclude_bands:
        exclude_keys.append(r)

infiles = sorted([x for x in os.listdir() if (x.endswith('10kb.txt') and x.startswith(keyword))])
window_dict = {r:0 for r in W}

def return_arr(infile,exclude_chroms=['chrX','chrY','chrM'],exclude_keys=exclude_keys,wd=window_dict):
    df = pd.read_csv(infile,sep='\t',usecols=['chrom','start','end','depth'])
    df = df[~df['chrom'].isin(exclude_chroms)]
    df['key'] = df.chrom + ':' + df.start.astype(str) + '-' + df.end.astype(str)
    df = df[~df['key'].isin(exclude_keys)]
    df = df[['key','depth']]
    keys = df['key'].tolist()
    vals = df['depth'].tolist()
    dict_x = wd.copy()
    if sum(list(dict_x.values()))!=0:
        print('error!')
    for k in keys:
        if k not in dict_x:
            print(k)
    dict_1 = dict(zip(keys, vals))
    dict_x.update(dict_1)
    vals = list(dict_x.values())
    return array(vals)

arrs = []
for i,f in enumerate(infiles):
    arr = return_arr(f)
    print(len(arr))
    arrs.append(arr)
arr = array(arrs,dtype=float).T
print(arr.shape)
savez(outfile,arr=arr)
#df.to_csv('merged_table_without_correction.csv',index=False)
