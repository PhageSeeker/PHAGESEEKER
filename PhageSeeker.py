#!/usr/bin/env python
# coding: utf-8


#libary

import os
import pandas as pd
import argparse
import subprocess


#argparse

parser = argparse.ArgumentParser()

parser.add_argument("-c","--crispr_path", help="Insert the whole directory path to reach the folder with the crispr detect file")
parser.add_argument("-p", "--phage_path", help="Insert the whole directory path to reach the folder with the Phage dataset file")
parser.add_argument("-o", "--folder_name", help="Insert the name of the folder where the spacer file will be stored")

args = parser.parse_args()
CRISPR_folder_path = args.crispr_path
PHAGES_input = args.phage_path
folder_name = args.folder_name


print('PhageSeeker at your service')
print('Code initialization....')



#extract spacer

def extract_spacer(path):
    
    col = ['contg','soft','spaer', 'start','end','length','direction','dot','ID']
    
    data = pd.read_csv(path, sep ='\t', names=col)
    
    id_list = []

    note_list = []


    for k in range(len(data)):
        split_id = data['ID'][k].split(';')
    
        id_list.append(list(filter(lambda word: word.startswith('ID'), split_id)))
    
        note_list.append(list(filter(lambda word: word.startswith('Note'), split_id)))
    

    id_list = [i[0] for i in id_list]
    
    note_list = [i[0][5:] for i in note_list]
    
    SR_list = []

    for i in id_list:
        a = i.split('_')
        SR_list.append(a[1])


    data['DNA'] = note_list

    data['S/R'] = SR_list

    new_data_list = ['contg','S/R','start','end','length','DNA']

    data = data[new_data_list]

    spacer = []


    for i in data.index:
    
        if 'SPACER' in data['S/R'][i]:
        
            spacer.append(i)
            

    
    filt_data = data.iloc[spacer,:]
    
    
    return filt_data


#create_spacer_file

def create_spacer_file(path,file_name,data): 
    
    path_f = os.path.join(folder_path,file_name)
    
    
    with open(path_f, 'w') as file:   
        
        for i in data.index:
            
            sequence = filt_data['DNA'][i]
            
            name = filt_data['contg'][i]
            
            fasta_entry = f'>{name}\n{sequence}\n'
            file.write(fasta_entry)
            
    return fasta_entry



#run_bash_command

def run_bash_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode('utf-8'), stderr.decode('utf-8'), process.returncode



#blast_and_plot

def blast_and_plot(PHAGES_input,file_name,path):

    command = "blastn -query "+ path + folder_name +'/' + file_name +" -subject "+ PHAGES_input +" -outfmt 6 -out CRISPRresult.txt"
    
    run_bash_command(command)

    names = ['MAG','phages ID', 'matching_%','spacer_length', 'mismatch',
             'gap_open','q_start','q_end','s_start','s_end','evalue','bitscore']
    
    
    result = pd.read_csv('/home/pol/CRISPRresult.txt', sep ='\t', names=names)
    
    return result


#iterate all the function 

for subdir, dirs, files in os.walk(CRISPR_folder_path):
    
    folder_path = str(subdir) + folder_name
    
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        
    tot_result = pd.DataFrame() 
    
    for file in files:
        
        path_t = os.path.join(subdir, file)
        
        
        filt_data = extract_spacer(path_t) #
        
        print('done filtering')
        
        
        file_name = 'spacer' + str(file)
        
        create_spacer_file(subdir,file_name,filt_data) #
        
        print('stored the spacer')
        
        result = blast_and_plot(PHAGES_input,file_name,subdir) #
        
        print('made the comparison')
        
        if len(result) > 0:
            
            tot_result = pd.concat([tot_result, result])
            
        
file_path = CRISPR_folder_path + folder_name + '/PhageDetect.csv'
tot_result.to_csv(file_path, index=False)       
        
print('Finish')

print('you can find your result in ' + folder_name)
        






