#!/usr/bin/env python
# coding: utf-8

import os
import glob
import pandas as pd
import warnings
warnings.simplefilter('ignore', FutureWarning)
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def tRNAcorrelation(input1, input2, output_folder):
    print('Input data1: '+input1)
    print('Input data2: '+input2)
    print('Output folder: '+output_folder+'\n')

    dirname = output_folder+'/'
    os.makedirs(dirname, exist_ok=True)

    #Extract file name. 
    filepath1 = input1+'/*profile.txt'
    filename1 = glob.glob(filepath1)

    filepath2 = input2+'/*profile.txt'
    filename2 = glob.glob(filepath2)

    #sort the list
    filename1.sort()
    filename2.sort()

    print('detected files in input 1: \n'+str(filename1)+'\n')
    print('detected files in input 2: \n'+str(filename2)+'\n')

    compiled_data=pd.DataFrame()
    
    #read reactivity data
    for in1, in2 in zip(filename1, filename2):
        data1 = pd.read_csv(in1, sep='\t', header='infer')
        reactivity1 = data1[['Reactivity_profile']]
        data2 = pd.read_csv(in2, sep='\t', header='infer')
        reactivity2 = data2[['Reactivity_profile']] 
        compiled_reactivity=pd.concat([reactivity1, reactivity2], axis=1)
        compiled_data=compiled_data.append(compiled_reactivity)
    df = compiled_data.set_axis(['input1', 'input2'], axis='columns')
    df_path=dirname+'A_correlation.csv'
    df = df.reset_index()
    del df['index']
    df.to_csv(df_path, index=None, header=True)
    print(df)
    
    #fig
    fig_name=dirname+'A_correlation_plot.png'
    fig=sns.jointplot(x='input1', y='input2', data=df, xlim=(-0.2, 1), ylim=(-0.2, 1), kind='reg',joint_kws={'line_kws':{'color':'grey'}})
    fig.set_axis_labels('Mutation rate (input1)', 'Mutation rate (input2)')
    plt.savefig(fig_name, dpi=300, bbox_inches='tight')

    #plt.clf() 
   


def main():
    parser=argparse.ArgumentParser(description='')
    parser.add_argument('-i1', '--input1', type=str, required=True, help='Specify input 1.')
    parser.add_argument('-i2', '--input2', type=str, required=True, help='Specify input 2.')
    parser.add_argument('-o', '--output_folder', type=str, default='Output', help='Results are provided in this folder. Specify a folder name.')

    args=parser.parse_args()

    tRNAcorrelation(args.input1, args.input2, args.output_folder)

if __name__ == '__main__':
    main()