#!/usr/bin/env python3
# coding: utf-8


'''
GDL Code for cointegrates

2024 03 30

Concatenates all GZ files into one file named for the barcode, and moves it up one tier in folder architecture (out of the barcode sub folder)

'''

import os
from Cointegrate_functions import *

directory = '/Users/glampe/Desktop/Columbia/Sternberg_Data/NGS/20240604_Cointegrates' # directory with csv input and folder named "reads" with all the barcode folders inside

output_folder = "./Reads"
reads_folder_path = "./Reads"

csv_input = '20240604_cointegrate_check.csv' # CSV input file that has barcode number, sample description, target sequence, TnSeq1, TnSeq2, and donor backbone for cointegrate checks (depending on primer, TnSeqs may vary)

'''

File format should be: 
Barcode, sample_description, target_sequence, tn_seq_1, tn_seq_2, cointegrate_seq

'''

output_csv = '20240604_output.csv'

QScore = 8 # QScore threshold for read filtering

kmer = 20 #kmer size to screen for target sequences, 20 is typically what I do

hdist = 2 # hamming distance for kmer search

# keep track of file paths to make sure it's all gonna path okay

os.chdir(directory)

for item in ['Read_histograms/','Cointegrate_outputs/','Cointegrate_outputs/QScore_filter/','Cointegrate_outputs/Insertions/','Cointegrate_outputs/Cointegrates/']:
    if not os.path.isdir(item): os.mkdir(item)

sample_dict = {}

process_input_file(csv_input,sample_dict)

concatentate_fastqgz(reads_folder_path,output_folder)

os.system('cls' if os.name == 'nt' else 'clear')    # great to clear terminal if you need to

for sample in sample_dict: analyze(sample,sample_dict,QScore,directory,kmer,hdist)

write_output(sample_dict,output_csv)





'''

CO INTEGRATES SHOULD BE DISCOUNTED FROM OTHER PCR COUNT RIGHT?!?! I MEAN IF EVERY SINGLE INSERTION IS COINTEGRATE, YOU WILL GET .5 SINCE YOU HAVE THE OTHER PCR

'''
