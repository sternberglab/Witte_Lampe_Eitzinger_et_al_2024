#!/usr/bin/env python3
# coding: utf-8
### Code by George Lampe, adapted from Sanne Klompe and Chris Acree ###
### requires input from Chris's input.csv format, ALL SEQUENCES MUST BE IN UPPER CASE ###
### updated 2022/10/12 ###

from Bio import SeqIO
import fnmatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from Bio import motifs
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import re
from pathlib import Path
from itertools import islice
import shutil
import os
import csv
import time
from datetime import datetime
import xlsxwriter

# set directory to wherever your reads are
os.chdir('/Users/glampe/Desktop/Columbia/Sternberg_Data/NGS/20231017_Isaac_RC_Test/Reads/')

#set output name at excel_out_name
excel_out_name = "/Users/glampe/Desktop/Columbia/Sternberg_Data/NGS/20231017_Isaac_RC_Test/test.xlsx"

# put path to input file, make sure it's in the correct format
csv_input = '/Users/glampe/Desktop/Columbia/Sternberg_Data/NGS/20231017_Isaac_RC_Test/isaac_input.csv'

seq_dated = input("Do sequence files have the date appended? y/n  ")
if seq_dated == "y": sequencing_date = input("\nPlease enter sequencing date on file names in YYYYMMDD format: ")

if not os.path.isdir('Distributions/'):
    os.mkdir('Distributions/')

def SequenceTally(filename, line, row, codename, sample_name):  # OLD VERSION DOESN'T USE BBDUK main command function for filtering
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time = {}".format(current_time))
    start_time = time.time()
    print("Analyzing Sample {}".format(row[0]))
    sample_distribution = {

    }
    Tn_End = row[2][0:10]
    Tn_RC = Seq(Tn_End.upper()).reverse_complement()
    Description = row[1]
    genome = Seq(row[5].upper()).reverse_complement()
    crRNA = Seq(row[4].upper()).reverse_complement()
    Genomic_crRNA = genome.find(crRNA)
    Genomic_site = genome[Genomic_crRNA-30:Genomic_crRNA-20]
    Unintegrated_site = genome[Genomic_crRNA-60:Genomic_crRNA-50]

    total = 0
    Unintegrated_tally = 0
    Integrated_tally = 0
    Neither = 0
    Waste_read = 0

    for record in SeqIO.parse(filename, "fastq"):
        total += 1
        genome_read = record.seq.find(Genomic_site)
        insertion_read = record.seq.find(Tn_RC)
        unintegrated_read = record.seq.find(Unintegrated_site)
        if genome_read > -1:
            if insertion_read > -1:
                Integrated_tally += 1
                integration_distance = genome_read - insertion_read + 20
                
                if integration_distance not in sample_distribution:
                    sample_distribution[integration_distance] = 1
                else:
                    sample_distribution[integration_distance] += 1
            elif insertion_read == -1:
                if unintegrated_read > -1:
                    Unintegrated_tally += 1
                elif unintegrated_read == -1:
                    Neither += 1
        if genome_read == -1:
            Waste_read += 1
    logsheet.write(line, 0, Description)
    logsheet.write(line, 1, total)
    logsheet.write(line, 2, Unintegrated_tally)
    logsheet.write(line, 3, Integrated_tally)
    logsheet.write(line, 4, Neither)
    logsheet.write(line, 5, Waste_read)
    fig = plt.figure(figsize=(5, 5), dpi=120)
    plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_title('{}\nSpacer length: {}nt'.format(Description,len(crRNA)))
    ax1.set_xlabel('Distance (bp)')
    ax1.set_ylabel('Reads')
    ax1.set_xlim([43,57])
    ax1.set_xticks([45,50,55])
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    plt.bar(sample_distribution.keys(), sample_distribution.values(), color='#b45c4b', width=1)
    plt.savefig("Distributions/{}.svg".format(codename))
    plt.close(fig)
    print("Sample {} Completed \n".format(filename))
    print("Analyzed in {} seconds".format(time.time() - start_time))
    return

# def read_len(input_file):
#     with open(input_file, 'r') as output:
#         lines = output.readlines()
#         return(int(len(lines) / 4))

# def bbducky(input_file,output_file,sequence):
#     bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'outm=' + output_file, 'k=25', 'hdist=3', 'literal=' + sequence, 'rcomp=f', '-Xmx6g'], stderr=subprocess.PIPE, text=True)


# def SequenceTally(filename, line, row, codename, sample_name):  # main command function for filtering
#     now = datetime.now()
#     current_time = now.strftime("%H:%M:%S")
#     print("Current Time = {}".format(current_time))
#     start_time = time.time()
#     print("Analyzing Sample {}".format(row[0]))

#     # Tn_End = row[2][0:15]
#     Tn_End = row[6]
#     Description = row[1]
#     genome = Seq(row[5])
#     crRNA = Seq(row[4])
#     Genomic_crRNA = genome.find(crRNA)
#     # Genomic_site = genome[Genomic_crRNA+65:Genomic_crRNA+75]
#     Unintegrated_site = row[7]
#     # Unintegrated_site = genome[Genomic_crRNA+83:Genomic_crRNA+98]

#     filez = ['unintegrated.fastq','integrated.fastq']

#     # genomic_filtered_file = filez[0]
#     unintegrated_file = filez[0]
#     integrated = filez[1]

#     # bbducky(filename,genomic_filtered_file,str(Genomic_site))
#     bbducky(filename,unintegrated_file,Unintegrated_site)
#     bbducky(filename,integrated,Tn_End)

#     total = read_len(filename)

#     # Genomic_sites = read_len(genomic_filtered_file)

#     Unintegrated_tally = read_len(unintegrated_file)

#     Integrated_tally = read_len(integrated)

#     Genomic_sites = Unintegrated_tally + Integrated_tally

#     Neither = total - Unintegrated_tally - Integrated_tally

#     Waste_read = total - Genomic_sites



#     logsheet.write(line, 0, Description)
#     logsheet.write(line, 1, total)
#     logsheet.write(line, 2, Unintegrated_tally)
#     logsheet.write(line, 3, Integrated_tally)
#     logsheet.write(line, 4, Neither)
#     logsheet.write(line, 5, Waste_read)

#     for item in filez:
#         os.remove(item)

#     print("Sample {} Completed \n".format(filename))
#     print("Analyzed in {} seconds".format(time.time() - start_time))
#     return

excel_out_path = Path(excel_out_name)
if not excel_out_path.is_file():
    log = xlsxwriter.Workbook(excel_out_name)
    bold = log.add_format({'bold': True})
    logsheet = log.add_worksheet("SequenceTally Log")
    logsheet.write(0, 0, "Sample", bold)
    logsheet.write(0, 1, "Total Reads", bold)
    logsheet.write(0, 2, "Unintegrated_tally", bold)
    logsheet.write(0, 3, "Integrated_tally", bold)
    logsheet.write(0, 4, "Neither", bold)
    logsheet.write(0, 5, "Waste_reads", bold)


with open(csv_input, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    next(reader, None)  # skips header row
    line = 1
    for row in reader:
        sample_name = row[0]
        if seq_dated == 'y': codename = sequencing_date + "-" + row[0]
        else: codename = row[0]
        #print(row[2][0:10])
        # search for file in cwd based on Sample ID ('code')
        filename = 'none'
        for i in os.listdir('.'):
            if fnmatch.fnmatch(i, "{}.fastq".format(codename)):
                filename = i
                break
        if filename != 'none':
            SequenceTally(filename, line, row, codename, sample_name)
        else:
            print("WARNING - File Not Found For {}".format(codename))
        line += 1
    print("End csv")
log.close()
