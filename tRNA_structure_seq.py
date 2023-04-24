#!/usr/bin/env python
# coding: utf-8
#Copy and paste the ShapeMapper output folder and put it into a working directory.

import os
import glob
import numpy as np
import subprocess
from subprocess import PIPE
import shutil
import fnmatch
import pandas as pd
import math
import warnings
warnings.simplefilter('ignore', FutureWarning)
import argparse
import matplotlib.pyplot as plt

def tRNAstructureseq(shapemapper_data, output_folder):
    print('Input data: '+shapemapper_data)
    print('Output folder: '+output_folder)

    proc=subprocess.run("conda activate tRNAstructureseq", shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print(proc.stdout)
    print(proc.stderr)

    dirname = 'Structure_seq_analysis/'+output_folder+'/'
    os.makedirs(dirname, exist_ok=True)

    #Extract file name. 
    filepath = shapemapper_data+"/*profile.txt"
    filename = glob.glob(filepath)

    #sort the list
    filename.sort()

    #Remove GU reactivity.
    print('Removing GU reactivities...\n')

    for i in filename:
        data = pd.read_csv(i, sep='\t', header='infer')
        if 'Norm_profile' in data:
            file_name = i.replace(".txt", "_GU_rm.shape")
            print(file_name)
            selected_columns = data[["Nucleotide", "Sequence", "Norm_profile"]]
            GU_rm_data = selected_columns.copy()
            GU_rm_data['Mask'] = np.where((GU_rm_data['Sequence'] == 'G') | (GU_rm_data['Sequence'] == 'U'), '-999', GU_rm_data['Norm_profile'])
            GU_rm_data['Mask'] = GU_rm_data['Mask'].replace('nan', -999)
            print(GU_rm_data['Mask'])
            GU_rm_data.to_csv(file_name, columns=['Nucleotide', 'Mask'], header=False, index=False, sep='\t')
        else:
            file_name = i.replace(".txt", "_GU_rm.shape")
            print(file_name)
            selected_columns = data[["Nucleotide"]]
            GU_rm_data = selected_columns.copy()
            GU_rm_data['Mask'] = np.nan
            GU_rm_data['Mask'] = GU_rm_data['Mask'].fillna(-999)
            GU_rm_data.to_csv(file_name, columns=['Nucleotide', 'Mask'], header=False, index=False, sep='\t')

    #Obtain new file names
    filepath = shapemapper_data+"/*profile_GU_rm.shape"
    Reactivity_file = glob.glob(filepath)

    #sort the list
    Reactivity_file.sort()

    #Fasta name
    filepath_fasta = "Structure_seq_analysis/reference_fasta/*.fasta"
    fasta_file = glob.glob(filepath_fasta)
    fasta_file.sort()

    #RNA name
    RNA_name = fasta_file
    RNA_name = [fasta_file.replace(".fasta", "") for fasta_file in RNA_name]
    RNA_name = [fasta_file.replace("Structure_seq_analysis/reference_fasta/", dirname) for fasta_file in RNA_name]
    RNA_name.sort()

    #CT name
    filepath_CT = "Structure_seq_analysis/reference_ct/*accept.ct"
    CT_file = glob.glob(filepath_CT)
    CT_file.sort()

    #Structure prediction
    print('Predicting structures...\n')
    for shape_file, tRNA_name, fasta_name, CT_name in zip(Reactivity_file, RNA_name, fasta_file, CT_file):
        print("Started "+tRNA_name)
        argument1="Fold "+fasta_name+" "+tRNA_name+"-Fold.ct -dms "+shape_file+" -mfe"
        print("argument1: " +argument1)
        proc=subprocess.run(argument1, shell=True, stdout=PIPE, stderr=PIPE, text=True)
        print(proc.stdout)
        print(proc.stderr)
        argument2="ct2dot "+tRNA_name+"-Fold.ct 1 "+tRNA_name+".dot"
        print("argument2: "+argument2)
        proc=subprocess.run(argument2, shell=True, stdout=PIPE, stderr=PIPE, text=True)
        print(proc.stdout)
        print(proc.stderr)
        argument3="draw "+tRNA_name+"-Fold.ct "+tRNA_name+"-Fold.ps -n 1 -s "+shape_file+" -l"
        print("argument3: "+argument3)
        proc=subprocess.run(argument3, shell=True, stdout=PIPE, stderr=PIPE, text=True)
        print(proc.stdout)
        print(proc.stderr)
        argument4="scorer "+tRNA_name+"-Fold.ct "+CT_name+" "+tRNA_name+"-MFE-scorer.txt"
        print("argument4: "+argument4)
        proc=subprocess.run(argument4, shell=True, stdout=PIPE, stderr=PIPE, text=True)
        print(proc.stdout)
        print(proc.stderr)
        print("----------------done----------------")
        print()

    path_scorer = dirname+'*-MFE-scorer.txt'
    argument5='cat '+path_scorer+'> '+dirname+'A_output_MFE_summary.txt'
    proc=subprocess.run(argument5, shell=True, stdout=PIPE, stderr=PIPE, text=True)
    print(proc.stdout)
    print(proc.stderr)
    print('saving the summarized data as A_output_MFE_summary.txt\n\n')

    print('Analyzing prediction acuracy\n\n')
    path_MFE_summary=dirname+'A_output_MFE_summary.txt'
    df=pd.read_table(path_MFE_summary, encoding="utf_8", sep=":", header=None, names=["A","B"])

    #dataframe
    summarized_data=pd.DataFrame(columns=['tRNA','Sensitivity','PPV', 'MCC'])
    summarized_data

    #Read sample name, Sensitivity, PPV, and append the information to the CSV file
    for i in range(0,len(df),5):
        tRNA_name=df.iloc[i]["B"]
        tRNA_name=tRNA_name.replace("-accept.ct", "")
        tRNA_name=tRNA_name.replace("  ", "")
        tRNA_name=tRNA_name.replace('reference_ct/', '')
        Sensitivity=df.iloc[i+3]["B"]
        Sensitivity=Sensitivity[-7:-1]
        #print("Sensitivity: "+Sensitivity)
        PPV=df.iloc[i+4]["B"]
        PPV=PPV[-7:-1]
        #print("PPV: "+PPV)
        MCC=math.sqrt(float(Sensitivity)*float(PPV))
        #print("MCC: "+str(MCC))
        summarized_data=summarized_data.append({"tRNA":tRNA_name, "Sensitivity":Sensitivity,
                                            "PPV":PPV, "MCC":float(MCC)}, ignore_index=True)
    print('Structure prediction accuracy:')
    print(summarized_data)
    print('saving dataframe as a csv file (A_Summarized_MFE_data.csv)')

    #Output dataframe as a csvfile
    path_summarized_data=dirname+'A_Summarized_MFE_data.csv'
    summarized_data.to_csv(path_summarized_data, index=False, header=True)

    #draw fig
    fig_name=dirname+'A_tRNAstructure_fig.png'
    print('Drawing figure...\n')
    fig_size = len(summarized_data['tRNA'])/5
    x=summarized_data['tRNA']
    y=summarized_data['MCC']
    fig=plt.figure() 
    fig.set_figwidth(fig_size)
    plt.bar(x, y, color='#007ACC', alpha=0.2, linewidth=5)
    plt.plot(x, y, "o", markersize=5, color='#007ACC', alpha=0.6)
    plt.ylabel('Prediction accuracy', fontsize=12, color = '#333F4B')
    plt.xlabel('')
    plt.xticks(rotation=90)
    plt.savefig(fig_name, dpi=300, bbox_inches='tight')
    plt.clf() 
    print('-------------------Done!-------------------\n')      

def main():
    parser=argparse.ArgumentParser(description='This script performs tRNA Structure-seq analysis using shapemapper outputs.')
    parser.add_argument('-s', '--shapemapper', type=str, required=True, help='Specify a name of shapemapper output.')
    parser.add_argument('-o', '--output_folder', type=str, default='Output_structure_prediction', help='Results are provided in this folder. Specify a folder name.')

    args=parser.parse_args()

    tRNAstructureseq(args.shapemapper, args.output_folder)

if __name__ == '__main__':
    main()