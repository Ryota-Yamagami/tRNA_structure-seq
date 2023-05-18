#!/usr/bin/env python
# coding: utf-8

import os
import time
import glob
import subprocess
from subprocess import PIPE
import shutil
import fnmatch
import pandas as pd
import warnings
warnings.simplefilter('ignore', FutureWarning)

def gtRNA_setup():
    #### generate output folder ####
    out1='Structure_seq_analysis/reference_fasta'
    os.makedirs(out1, exist_ok=True)
    out2='Structure_seq_analysis/reference_dot'
    os.makedirs(out2, exist_ok=True)
    out3='Structure_seq_analysis/reference_ct'
    os.makedirs(out3, exist_ok=True) 
    out4='MaP_analysis/reference_fasta_ShapeMapper'
    os.makedirs(out4, exist_ok=True) 

    #### read ss.sort and name_map.txt files
    file_names=glob.glob('*-tRNAs/*')
    
    if any(var.endswith('-confidence-set.ss') for var in file_names):
        file_ss=fnmatch.filter(file_names, '*confidence-set.ss')
        file_sort=file_ss[0]+'.sort'
        shutil.copyfile(file_ss[0], file_sort)
    else:
        pass

    ss_filepath='*-tRNAs/*.sort'
    ss_filename=glob.glob(ss_filepath)
    print(ss_filename)
    ss_df=pd.read_csv(ss_filename[0], sep='\s+', usecols=[0, 1, 2], names=['tRNA_index', 'a', 'b'])
    print(ss_df)


    name_filepath='*-tRNAs/*map.txt'
    name_filename=glob.glob(name_filepath)
    name_df=pd.read_csv(name_filename[0], sep='\t', header='infer')
    SE_id=list(name_df['tRNAscan-SE_id'])
    tRNA_name=list(name_df['GtRNAdb_id'])

    for var, tRNA in zip(SE_id, tRNA_name):
        if 'NNN' in tRNA:
            pass

        elif tRNA[-2:] == '-1':
            print('Searching ', var)
            id_index=ss_df.loc[ss_df['tRNA_index'].str.fullmatch(var)]

            if id_index.empty:
                print("did't find ", var, " : skipping this unconfident tRNA.\n")
                pass 
            
            else:                
                print('Found at index:', id_index['tRNA_index'], '\n')
                index_num=id_index.index.values
                index_seq=index_num + 4
                index_ss=index_num + 5
                discrim_num=index_num + 2
                discrim=ss_df.iat[discrim_num[0], 1]

                if 'intron:' in discrim:
                    index_seq=index_num + 5
                    index_ss=index_num + 6
                    print('tRNA name: ', tRNA)
                    seq=ss_df.iat[index_seq[0], 1]
                    ss=ss_df.iat[index_ss[0], 1]
                    ss=ss.replace('>', '(').replace('<', ')')
                    print(seq)
                    print(ss, '\n')
                    print('###generating fasta and dot bracket files###')
                    fa_file_name=out1+'/'+tRNA[5:]+'.fasta'
                    fa_seq_name='>'+tRNA[5:]
                    dot_file_name=out2+'/'+tRNA[5:]+'-accept.dot'
                    
                    with open(fa_file_name, 'w') as f_fa:
                        f_fa.write(fa_seq_name+'\n')
                        f_fa.write(seq+'\n')
                    print('saved as a fasta format')

                    with open(dot_file_name, 'w') as f_dot:
                        f_dot.write(fa_seq_name+'\n')
                        f_dot.write(seq+'\n')
                        f_dot.write(ss+'\n')
                    print('saved as a dot bracket format')
                    
                    ct_file_name=out3+'/'+tRNA[5:]+'-accept.ct'
                    argument="dot2ct "+dot_file_name+' '+ct_file_name
                    proc=subprocess.run(argument, shell=True, stdout=PIPE, stderr=PIPE, text=True)
                    print('saved as a ct format\n###done!###\n')
                
                else:
                    print('tRNA name: ', tRNA)
                    seq=ss_df.iat[index_seq[0], 1]
                    ss=ss_df.iat[index_ss[0], 1]
                    ss=ss.replace('>', '(').replace('<', ')')
                    print(seq)
                    print(ss, '\n')
                    print('###generating fasta and dot bracket files###')
                    fa_file_name=out1+'/'+tRNA[5:]+'.fasta'
                    fa_seq_name='>'+tRNA[5:]
                    dot_file_name=out2+'/'+tRNA[5:]+'-accept.dot'
                    
                    with open(fa_file_name, 'w') as f_fa:
                        f_fa.write(fa_seq_name+'\n')
                        f_fa.write(seq+'\n')
                    print('saved as a fasta format')

                    with open(dot_file_name, 'w') as f_dot:
                        f_dot.write(fa_seq_name+'\n')
                        f_dot.write(seq+'\n')
                        f_dot.write(ss+'\n')
                    print('saved as a dot bracket format')
                    
                    ct_file_name=out3+'/'+tRNA[5:]+'-accept.ct'
                    argument="dot2ct "+dot_file_name+' '+ct_file_name
                    proc=subprocess.run(argument, shell=True, stdout=PIPE, stderr=PIPE, text=True)
                    print('saved as a ct format\n###done!###\n')
            
        
        else:
            pass
    
    #Merge fasta file
    argument="cat Structure_seq_analysis/reference_fasta/* > MaP_analysis/reference_fasta_ShapeMapper/reference_shapemapper.fa"
    proc=subprocess.run(argument, shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print('saved reference fasta for shapemapper \n###done!###\n')



def main():
    gtRNA_setup()

if __name__ == '__main__':
    t1 = time.time()
    main()
    t2 = time.time()
    elapsed_time = t2-t1
    print(f"Elapsed_time: {elapsed_time} sec")