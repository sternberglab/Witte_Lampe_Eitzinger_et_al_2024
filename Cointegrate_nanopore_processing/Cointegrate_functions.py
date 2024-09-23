#!/usr/bin/env python3
# coding: utf-8


'''
GDL Code for cointegrates

2024 03 30

Concatenates all GZ files into one file named for the barcode, and moves it up one tier in folder architecture (out of the barcode sub folder)

'''

import csv
import glob
import subprocess
import os
import shutil
import gzip

# keep track of file paths to make sure it's all gonna path okay

def bbduk_func(input, output, seq, kmer, hdist,output_unmatched = [False]):
    if not output_unmatched[0]:
    	return subprocess.run(['bbduk.sh', 'in1={}'.format(input), 'outm={}'.format(output), 
        	'k={}'.format(kmer), 'hdist={}'.format(hdist), 'literal={}'.format(seq), '-Xmx6g'], stderr=subprocess.PIPE, text=True)
    else:
    		return subprocess.run(['bbduk.sh', 'in1={}'.format(input), 'out={}'.format(output_unmatched[1]), 'outm={}'.format(output), 
    	    	'k={}'.format(kmer), 'hdist={}'.format(hdist), 'literal={}'.format(seq), '-Xmx6g'], stderr=subprocess.PIPE, text=True)

def process_input_file(csv_input,sample_dict): # processes input CSV file to get all samples to anlayze and their requisite sequences/descriptions
    with open(csv_input, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader) # Skip header row

        for row in reader:
            sample_dict[row[0]]= {
                'barcode': row[0],
                'name': row[1],
                'target_seq': row[2],
                'seq1': row[3],
                'seq2': row[4],
                'cointegrate': row[5]
            }

def concatentate_fastqgz(reads_folder_path,output_folder): #concatenates all gz files in each barcode folder within "reads" to then move one folder up and rename file to be the barcode

    # Iterate through all folders within the "Reads" folder
    for folder_name in os.listdir(reads_folder_path):
        folder_path = os.path.join(reads_folder_path, folder_name)
        if os.path.isdir(folder_path):
            # Find all FASTQ.gz files in the subfolder
            fastq_files = glob.glob(os.path.join(folder_path, "*.fastq.gz"))
            if fastq_files:
                # Concatenate all fastq.gz files
                cat_command = ["cat"] + fastq_files
                cat_process = subprocess.Popen(cat_command, stdout=subprocess.PIPE)
                
                # Write the concatenated output to a new file
                output_filename = os.path.join(output_folder, folder_name + ".fastq.gz")
                with open(output_filename, "wb") as output_file:
                    output_file.write(cat_process.stdout.read())
                
                print("Concatenation completed for folder:", folder_name)
                
                # Delete all files in the folder
                for file_name in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, file_name)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                print("Files deleted in folder:", folder_name)
                
                # Delete the parent folder
                shutil.rmtree(folder_path)
                print("Parent folder deleted:", folder_path)

def get_file_size(file_path,directory): # returns integer of number of reads in input gz file double check that this is a right method for nanopore fastQ files

	try:
	    with gzip.open(file_path, 'rt') as output:
	        read_count = output.readlines()
	    return int(len(read_count) // 4)  # Divide by 4 because each read consists of 4 lines
	except Exception as e:
	    print(f"An error occurred: {e}")
	    return 0

def qScore_filter(sample, qscore,sample_dict,directory): #filters for QScore
    input_path = 'Reads/{}.fastq.gz'.format(sample_dict[sample]['barcode'])
    output_path = 'Cointegrate_outputs/QScore_filter/{}_qscore.fastq.gz'.format(sample_dict[sample]['barcode'])

    bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_path, 'out=' + output_path, 'maq={}'.format(qscore),'lhist=Read_histograms/{}_lengths.txt'.format(sample_dict[sample]['barcode']),'-Xmx6g'], stderr=subprocess.PIPE, text=True)

    sample_dict[sample]['input_reads'] = get_file_size(input_path,directory)
    sample_dict[sample]['qscore_reads'] = get_file_size(output_path,directory)

    # with open('Cointegrate_outputs/QScore_filter/{}_QscoreFiltering_bbduk.txt'.format(sample_dict[sample]['barcode']),'w') as f: f.write(bbduk_process.stderr)

def filter_insertions(sample,sample_dict,kmer,hdist,directory): #filters for reads that have insertion events'
    input_path = 'Cointegrate_outputs/QScore_filter/{}_qscore.fastq.gz'.format(sample_dict[sample]['barcode'])
    target_filter = 'Cointegrate_outputs/Insertions/{}_targetseq.fastq.gz'.format(sample_dict[sample]['barcode'])
    seq1_filter = 'Cointegrate_outputs/Insertions/{}_seq1.fastq.gz'.format(sample_dict[sample]['barcode'])
    output_path = 'Cointegrate_outputs/Insertions/{}_insertions.fastq.gz'.format(sample_dict[sample]['barcode'])

    bbduk_process_1 = bbduk_func(input_path, target_filter, sample_dict[sample]['target_seq'], kmer, hdist)
    bbduk_process_2 = bbduk_func(target_filter, seq1_filter, sample_dict[sample]['seq1'], kmer, hdist)
    bbduk_process_3 = bbduk_func(seq1_filter, output_path, sample_dict[sample]['seq2'], kmer, hdist)

    for file in [target_filter,seq1_filter]: os.remove(file)
    sample_dict[sample]['insertions']=get_file_size(output_path,directory)

    # with open('Cointegrate_outputs/Insertions/{}_insertion_filterings1_bbduk.txt'.format(sample_dict[sample]['barcode']),'w') as f: f.write(bbduk_process_1.stderr)
    # with open('Cointegrate_outputs/Insertions/{}_insertion_filterings2_bbduk.txt'.format(sample_dict[sample]['barcode']),'w') as f: f.write(bbduk_process_2.stderr)
    # with open('Cointegrate_outputs/Insertions/{}_insertion_filterings3_bbduk.txt'.format(sample_dict[sample]['barcode']),'w') as f: f.write(bbduk_process_3.stderr)

def filter_cointegrates(sample,sample_dict,kmer,hdist,directory): # filters for insertion events that additionally have donor backbone (indicating cointegrate)
    input_path = 'Cointegrate_outputs/Insertions/{}_insertions.fastq.gz'.format(sample_dict[sample]['barcode'])
    output_path = 'Cointegrate_outputs/Cointegrates/{}_cointegrates.fastq.gz'.format(sample_dict[sample]['barcode'])
    output_unmatched_path = 'Cointegrate_outputs/Insertions/{}_simple_insertions.fastq.gz'.format(sample_dict[sample]['barcode'])

    bbduk_process = bbduk_func(input_path, output_path, sample_dict[sample]['cointegrate'], kmer, hdist,output_unmatched = [True,output_unmatched_path])

    sample_dict[sample]['cointegrate_reads'] = get_file_size(output_path,directory)
    sample_dict[sample]['simple_insertion_reads'] = get_file_size(output_unmatched_path,directory) - sample_dict[sample]['cointegrate_reads']


    # with open('Cointegrate_outputs/Cointegrates/{}_cointegrates_bbduk.txt'.format(sample_dict[sample]['barcode']),'w') as f: f.write(bbduk_process.stderr)

def write_output(sample_dict,output_csv):

	with open(output_csv, 'w', newline='') as stats_file: # Open stats file for writing
	    writer = csv.writer(stats_file)
	    
	    writer.writerow(['Barcode', 'Sample','Target sequence', 'Sequence 1', 'Sequence 2', 'Cointegrate Sequence',
	    'Input Reads', 'Quality filter reads','Total insertion reads','Simple Insertions','Cointegrates']) # Write header row

	    for sample in sample_dict:
	    	writer.writerow(
	    	    [sample_dict[sample]['barcode'], sample_dict[sample]['name'], sample_dict[sample]['target_seq'], sample_dict[sample]['seq1'], sample_dict[sample]['seq2'], sample_dict[sample]['cointegrate'],
	    	    sample_dict[sample]['input_reads'], sample_dict[sample]['qscore_reads'], sample_dict[sample]['insertions'], sample_dict[sample]['simple_insertion_reads'], sample_dict[sample]['cointegrate_reads']])

def analyze(sample,sample_dict,QScore,directory,kmer,hdist):
    qScore_filter(sample,QScore,sample_dict,directory)
    filter_insertions(sample,sample_dict,kmer,hdist,directory)
    filter_cointegrates(sample,sample_dict,kmer,hdist,directory)
    print('\n{} analyzed\n'.format(sample))

