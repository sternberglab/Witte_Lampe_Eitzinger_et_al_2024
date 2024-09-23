#!/usr/bin/env python3
# coding: utf-8

'''
Written by GDL 2023/03/20

Processes reads based on containing correct transposon end, removing pDonor contaminating reads and control spike-in reads, as well as self-targeting reads into the CRISPR array

Processing done by BBDuk

If crRNA included, sorts on/off-target integration events by filtering based on reads aligning to 100bp window downstream of crRNA

Filters mapped/unmapped, on/off-target, as well as crRNA mapped reads into separate FASTA files

Plots genome-wide alignment data, zooms in for bottom x% aligned data, and plots T-RL vs T-LR orientation (when applicable)

MUST PERFORM UMI PROCESSING BEFOREHAND IF UMIs WERE USED IN THE TAGMENTATION PROCESS

READ THE PARAMS FILE! THINGS TO CUSTOMIZE THERE

'''

import time
from datetime import datetime

import csv
import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import glob

import shutil

from parameters import *
from pipeline.sample_processing import *
from pipeline.alignment import *


from pipeline.plotting import *
# from pipeline.plotting_log_TnpA import *

from pipeline.normalize import *

start_time = time.time()

os.chdir(params['Directory'])
code_path = "{}/code/".format(params['Prefix'])

for item in [params['Prefix'],intermediate_path,code_path]:
    if not os.path.isdir(item): os.mkdir(item)

for item in ["TagTn.py","parameters.py"]: os.system('cp {}/{} {}/code/{}'.format(params['Code_path'],item, params['Prefix'],item))
os.system('cp {}/pipeline {}/code/pipeline -r'.format(params['Code_path'],params['Prefix']))    

if params['Run_processing']:
    process_time = time.time()
    with open(params['info_file'], 'r') as file_list:
        reader = csv.reader(file_list)
        next(reader) # Skip header row
        with open('{}/{}_Processing_statistics.csv'.format(params['Prefix'],params['Prefix']), 'w', newline='') as stats_file: # Open stats file for writing
            writer = csv.writer(stats_file)
            
            writer.writerow(['Sample', 'Description','System', 'Genome', 'crRNA',
            'Input Reads', 'Quality filter reads','Percent Yield','Tn-End containing reads','Percent Yield', 
            'pDonor reads', 'Percent pDonor',
            'Spike-in Reads','Percent Spike-in','Plasmid Mapping Reads',
            'Self-targeting Reads','Percent Self-targeting',
            'Reads for genome mapping','Weblogo','Minicircle Reads','Analyze off-targets']) # Write header row
            
            for row in reader:
                sample_metainfo = generate_sample_metainfo(row, params['sequencing_date'])
                os.system('mkdir {}{}/'.format(intermediate_path,sample_metainfo['Sample']))
                if len(sample_metainfo['off_targets'])>0: sample_metainfo['analyze_offs'] = 'Yes'
                else: sample_metainfo['analyze_offs'] = 'No'
                print("\n------\n\nProcessing {}\n\n------\n".format(sample_metainfo['Sample']))

                if len(sample_metainfo['crRNA'])>0: # writes crRNA file in alignment folder to map to genome
                    crRNA_fasta = SeqRecord(Seq(sample_metainfo['crRNA']), id = 'crRNA', description = sample_metainfo['Sample'])
                    if not os.path.exists('{}/alignment_input/crRNA_fasta_files/'.format(params['Prefix'])): os.makedirs('{}/alignment_input/crRNA_fasta_files/'.format(params['Prefix']))
                    SeqIO.write(crRNA_fasta, '{}/alignment_input/crRNA_fasta_files/{}.fasta'.format(params["Prefix"],sample_metainfo['Sample']), "fasta")
                
                if len(sample_metainfo['off_targets'])>0: # writes off-target crRNA file in alignment folder to map to genome
                    off_targets_fasta = []
                    for item in sample_metainfo['off_targets']:
                        off_target_site = SeqRecord(Seq(item), id = 'offtarget_{}'.format(item), description = sample_metainfo['Sample'])
                        off_targets_fasta.append(off_target_site)
                    if not os.path.exists('{}/alignment_input/offtarget_fasta_files/'.format(params['Prefix'])): os.makedirs('{}/alignment_input/offtarget_fasta_files/'.format(params['Prefix']))
                    SeqIO.write(off_targets_fasta, '{}/alignment_input/offtarget_fasta_files/{}.fasta'.format(params['Prefix'],sample_metainfo['Sample']), "fasta")
                
                filter_stats = quality_filter(sample_metainfo)
                if params['pDonor_UMIs']: UMI_stats = umi_tools(sample_metainfo)
                trimming_stats = Tn_finder_trimmer(sample_metainfo)
                pDonor_stats = pDonor_finder(sample_metainfo,trimming_stats)

                if sample_metainfo['System'] != "CAST" and sample_metainfo['System'] != "TnpA":

                    # simplest filtering, just looked for pDonor then nothing else

                    genome_reads = pDonor_stats[2]
                    spikein_stats = ['N/A','N/A']
                    CRISPR_stats = ["N/A", "N/A","N/A"]
                    plasmid_reads = 'N/A'
                
                elif sample_metainfo['System'] == "TnpA":

                    # For TnpA specifically, may need to tweak a bit

                    spikein_stats = spikein_finder(sample_metainfo,trimming_stats)
                    sample_metainfo['Minicircle_counts']=minicircle_finder(sample_metainfo,trimming_stats)
                    genome_reads = spikein_stats[2] - sample_metainfo["Minicircle_counts"][2]
                    plasmid_reads = 'N/A'
                    CRISPR_stats = ['N/A','N/A']

                elif sample_metainfo['System'] == "CAST":
                    spikein_stats = spikein_finder(sample_metainfo,trimming_stats)

                    if not params['self-target_mask'][0] and sample_metainfo['Genome']=='pace_ecoli': CRISPR_stats = ['NA','NA',spikein_stats[2]]
                    else: CRISPR_stats = self_targets(sample_metainfo,trimming_stats)
                    genome_reads = CRISPR_stats[2]
                    if sample_metainfo['Genome'] == 'human-noy' or sample_metainfo['Genome'] == 'human-noy_aavs1': 
                        plasmid_readout = plasmid_finder(sample_metainfo,trimming_stats)
                        genome_reads = plasmid_readout[1]
                        plasmid_reads = plasmid_readout[0]
                    else: plasmid_reads = "N/A"
                if not params['KEEP_INTERMEDIATES']: shutil.rmtree(intermediate_path+"{}/".format(sample_metainfo['Sample']))
                writer.writerow(
                    [sample_metainfo['Sample'], sample_metainfo['Description'], sample_metainfo['System'], 
                    sample_metainfo['Genome'], sample_metainfo['crRNA'],
                    filter_stats[0], filter_stats[1], filter_stats[2],
                    trimming_stats[1], trimming_stats[2], 
                    pDonor_stats[0], pDonor_stats[1], 
                    spikein_stats[0], spikein_stats[1], 
                    plasmid_reads,
                    CRISPR_stats[0], CRISPR_stats[1],
                    genome_reads,
                    sample_metainfo['Weblogo'],sample_metainfo['Minicircle_counts'][2],sample_metainfo['analyze_offs']])

    if not params['KEEP_INTERMEDIATES']: shutil.rmtree(intermediate_path)

    print("\nInitial processing complete\n\nProccessed in {} minutes".format((time.time() - process_time)//60))

elif not params['Run_processing']: print('\nRun_processing set to False, proceeding directly to alignment')

if not params['Run_alignment']: print('\nNot performing Alignment')

elif params['Run_alignment']:
    align_time = time.time()
    print("\nBeginning alignments and plotting\n\n\nChecking for proper index files\n")
    for item in index_genomes: 
        testing_index = check_index(alignment_algo, index_folder, index_genomes, item,cores_to_use)
        if not check_index: 
            print("Did not pass index check, exiting program")
            exit(1)
    
    print('\n\nIndexing complete \n\n')

    if not os.path.isdir('{}/{}_output/alignment_outputs'.format(params['Prefix'],alignment_algo)): os.makedirs('{}/{}_output/alignment_outputs'.format(params['Prefix'],alignment_algo))

    with open('{}/{}_Processing_statistics.csv'.format(params['Prefix'],params['Prefix']), 'r') as file_list:
        reader = csv.reader(file_list)
        next(reader)
        for row in reader: # Skip header row
            sample_metainfo = {
            'Sample': row[0],
            'Description': row[1],
            'System': row[2],
            'Genome': row[3].lower(),
            'crRNA': row[4],
            'off_targets': row[20],
            'Candidate alignment reads': int(row[17])
            }

            if sample_metainfo['Candidate alignment reads'] == '0':
                print("------\n------\n\nNo reads to analyze for {}!\n\n------\n------\n".format(sample_metainfo['Sample']))
                continue
            print("\n------\n\nAligning {} with {} cores\n\n------\n".format(sample_metainfo['Sample'],cores_to_use))

            index_path = generate_index_path(alignment_algo, index_folder, index_genomes, sample_metainfo['Genome'])

            align(index_path[1],sample_metainfo,cores_to_use,alignment_algo,params['is_paired_end'])
    print("\nAligned in {} minutes".format((time.time() - align_time)//60))


if not params['Run_plotting']: print("\n\nNot plotting alignments, program done")

if params['Run_plotting']:
    plot_time = time.time()
    max_spikein = 0
    max_norm_value = 0
    max_norm_values = []
    print("\n\nPerforming plotting")
    bins_dict = {}

    for item in ['','Distributions/','Raw/','Percent/','CSV_Files/','WeblogoFastas/','Weblogos/','Weblogos_Filtered/','Filtered_Raw/']:
        if not os.path.isdir("{}/Output_figures/{}".format(params['Prefix'],item)): os.mkdir("{}/Output_figures/{}".format(params['Prefix'],item))

    if params['TagTn_UMIs'] or params['pDonor_UMIs']:
        for item in ['','Distributions/','Raw/','Percent/','CSV_Files/','WeblogoFastas/','Weblogos/','Weblogos_Filtered/','Filtered_Raw/']:
            if not os.path.isdir("{}/UMI_output_figures/{}".format(params['Prefix'],item)): os.mkdir("{}/UMI_output_figures/{}".format(params['Prefix'],item))


    with open('{}/{}_Processing_statistics.csv'.format(params['Prefix'],params['Prefix']), 'r') as file_list:
        reader = csv.reader(file_list)
        next(reader)
        with open('{}/{}_Alignment_statistics.csv'.format(params['Prefix'],params['Prefix']), 'w', newline='') as stats_file: # Open stats file for writing
            writer = csv.writer(stats_file)
            writer.writerow(['Sample', 'Description','System', 'crRNA','Input Reads','Spike-in reads',
            'Mapped Reads', 'Unmapped Reads', 'Percent Yield',
            'On-target reads', 'Off-target reads', 'Percent On-target', 'Unique sites',
            'Percent T-RL',"Dual Reads","Misaligned Reads","T7 insertions","Self-targeting Reads","Mismatched Alignment","Filtered Mapped Reads","Filtered On","Filtered Off"])
            # NOTE: if the sample is not paired-end sequencing, misaligned won't mean anything!
            for row in reader: # Skip header row
                if int(row[5]) > max_spikein: max_spikein = int(row[5])
                sample_metainfo = {
                'Sample': row[0],
                'Description': row[1],
                'System': row[2],
                'Genome': row[3].lower(),
                'crRNA': row[4],
                'spike-in': row[12],
                'Candidate alignment reads': int(row[17]),
                'off_targets': row[20],
                'Weblogo':row[18]
                }
                if 'UMI' in sample_metainfo['Sample']: data = 'UMI'
                else: data = "read"
                if sample_metainfo['Candidate alignment reads'] == 0:
                    print("------\n------\n\nNo reads to analyze for {}!\n\n------\n------\n".format(sample_metainfo['Sample']))
                    writer.writerow(
                        [sample_metainfo['Sample'],sample_metainfo['Description'], sample_metainfo['System'], sample_metainfo['crRNA'],sample_metainfo['Candidate alignment reads'],
                        "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A","N/A","N/A"])
                    continue
                print("\n------\n\nPlotting {}\n\n------\n".format(sample_metainfo['Sample']))
                
                sample_metainfo.update(plotfig(sample_metainfo, TSD, alignment_algo,width,data,params['is_paired_end']))

                if sample_metainfo['Genome'] in ['pace_ecoli','ecoli']: print("\n----\n\nSample {} had {} insertions at the T7 locus\n\n----\n".format(sample_metainfo['Sample'],sample_metainfo['T7 insertions']))

                sample_metainfo['yield'] = str(100*sample_metainfo['mapped']/sample_metainfo['Candidate alignment reads'])
                if not sample_metainfo['ontarget']:
                    sample_metainfo['ontarget'] = "N/A"
                    sample_metainfo['offtarget'] = "N/A"
                    sample_metainfo['specificity'] = 'N/A'
                    sample_metainfo['RL'] = 'N/A'
                else:
                    sample_metainfo['specificity'] = str(100*sample_metainfo['ontarget']/(sample_metainfo['ontarget']+sample_metainfo['offtarget']))
                    sample_metainfo['specificity'] = sample_metainfo['specificity'].split(".")[0]+"."+sample_metainfo['specificity'].split(".")[1][:3]
                    sample_metainfo['RL'] = str(100*sample_metainfo['RL_counter']/(sample_metainfo['RL_counter']+sample_metainfo['LR_counter']))
                    sample_metainfo['RL'] = sample_metainfo['RL'].split(".")[0]+"."+sample_metainfo['RL'].split(".")[1][:3]
                writer.writerow(
                    [sample_metainfo['Sample'], sample_metainfo['Description'],sample_metainfo['System'], sample_metainfo['crRNA'],sample_metainfo['Candidate alignment reads'],sample_metainfo['spike-in'],
                    sample_metainfo['mapped'], (sample_metainfo['Candidate alignment reads'] - sample_metainfo['mapped']), 
                    sample_metainfo['yield'].split(".")[0]+"."+sample_metainfo['yield'].split(".")[1][:3],
                    sample_metainfo['ontarget'], sample_metainfo['offtarget'], sample_metainfo['specificity'], sample_metainfo['unique_sites'],
                    sample_metainfo['RL'],sample_metainfo['dual_reads'],sample_metainfo['misaligned'],sample_metainfo['T7 insertions'],sample_metainfo['self_target'],sample_metainfo['mismatched_alignment'],
                    sample_metainfo['filtered_total'],sample_metainfo['filtered_on'],sample_metainfo['filtered_off']])
                if sample_metainfo['Weblogo'] == 'Yes':
                    if sample_metainfo['unique_sites']!=0:
                        print("\n\nGenerating Fasta for Weblogo for Sample {}\n\n{} bases on each side are taken\n\n".format(sample_metainfo['Sample'],params['weblogo_flankseq']))
                        weblogo(sample_metainfo, alignment_algo)
                    else: print("\n\nNot generating weblogo for {}, no mapped integration sites\n\n".format(sample_metainfo['Sample']))
                max_norm_values.append(sample_metainfo['max_bin'])
                for item in max_norm_values:
                    if item*max_spikein>max_norm_value: max_norm_value = item*max_spikein
                bins_dict[sample_metainfo['Sample']]=sample_metainfo['plot_bins']

    if not params['spikein_norm']: print("\n\nNo spike-in normalization")
    else:
        normalize_time = time.time()
        print("\nBeginning Normalization\n")

        if not os.path.isdir(params['Prefix']+'/Output_figures/Normalized'): os.mkdir(params['Prefix']+"/Output_figures/Normalized")

        with open('{}/{}_Processing_statistics.csv'.format(params['Prefix'],params['Prefix']), 'r') as file_list:
            reader = csv.reader(file_list)
            next(reader)
            for row in reader:
                sample_metainfo = {
                'Sample': row[0],
                'Description': row[1],
                'System': row[2],
                'Genome': row[3].lower(),
                'crRNA': row[4],
                'spike-in': row[12],
                'Candidate alignment reads': int(row[17]),
                'off_targets': row[20],
                }

                if sample_metainfo['Candidate alignment reads'] == 0: continue
                else: sample_metainfo['plot_bins'] = bins_dict[sample_metainfo['Sample']]
                plotnormalizefig(sample_metainfo, sample_metainfo['plot_bins'],width,"normalized read",max_spikein,max_norm_value,alignment_algo)

        print("\nNormalized in {} minutes".format((time.time() - normalize_time)//60))

    print("\nPlotted in {} minutes".format((time.time() - plot_time)//60))

    if params['TagTn_UMIs'] or params['pDonor_UMIs']:
        print("\n\nPerforming UMI plotting")
        data = 'UMI'

        with open('{}/{}_Processing_statistics.csv'.format(params['Prefix'],params['Prefix']), 'r') as file_list:
            reader = csv.reader(file_list)
            next(reader)
            with open('{}/{}_UMI_Alignment_statistics.csv'.format(params['Prefix'],params['Prefix']), 'w', newline='') as stats_file: # Open stats file for writing
                writer = csv.writer(stats_file)
                writer.writerow(['Sample', 'Description','System', 'crRNA','Input Reads','Spike-in reads',
                'Mapped Reads', 'Unmapped Reads', 'Percent Yield',
                'On-target reads', 'Off-target reads', 'Percent On-target', 'Unique sites',
                'Percent T-RL',"Dual Reads","Misaligned Reads","T7 insertions","Self-targeting Reads","Mismatched Alignment","Filtered Mapped Reads","Filtered On","Filtered Off"])
                # NOTE: if the sample is not paired-end sequencing, misaligned won't mean anything!
                for row in reader: # Skip header row
                    if int(row[5]) > max_spikein: max_spikein = int(row[5])
                    sample_metainfo = {
                    'Sample': row[0],
                    'Description': row[1],
                    'System': row[2],
                    'Genome': row[3].lower(),
                    'crRNA': row[4],
                    'spike-in': row[12],
                    'Candidate alignment reads': int(row[17]),
                    'off_targets': row[20],
                    'Weblogo':row[18]
                    }
                    if sample_metainfo['Candidate alignment reads'] == 0:
                        print("------\n------\n\nNo reads to analyze for {}!\n\n------\n------\n".format(sample_metainfo['Sample']))
                        writer.writerow(
                            [sample_metainfo['Sample'],sample_metainfo['Description'], sample_metainfo['System'], sample_metainfo['crRNA'],sample_metainfo['Candidate alignment reads'],
                            "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A","N/A","N/A"])
                        continue
                    print("\n------\n\nPlotting {}\n\n------\n".format(sample_metainfo['Sample']))
                    
                    sample_metainfo.update(plotfig(sample_metainfo, TSD, alignment_algo,width,data,params['is_paired_end'],UMI_plotting=True))

                    if sample_metainfo['Genome'] in ['pace_ecoli','ecoli']: print("\n----\n\nSample {} had {} insertions at the T7 locus\n\n----\n".format(sample_metainfo['Sample'],sample_metainfo['T7 insertions']))

                    sample_metainfo['yield'] = str(100*sample_metainfo['mapped']/sample_metainfo['Candidate alignment reads'])
                    if not sample_metainfo['ontarget']:
                        sample_metainfo['ontarget'] = "N/A"
                        sample_metainfo['offtarget'] = "N/A"
                        sample_metainfo['specificity'] = 'N/A'
                        sample_metainfo['RL'] = 'N/A'
                    else:
                        sample_metainfo['specificity'] = str(100*sample_metainfo['ontarget']/(sample_metainfo['ontarget']+sample_metainfo['offtarget']))
                        sample_metainfo['specificity'] = sample_metainfo['specificity'].split(".")[0]+"."+sample_metainfo['specificity'].split(".")[1][:3]
                        sample_metainfo['RL'] = str(100*sample_metainfo['RL_counter']/(sample_metainfo['RL_counter']+sample_metainfo['LR_counter']))
                        sample_metainfo['RL'] = sample_metainfo['RL'].split(".")[0]+"."+sample_metainfo['RL'].split(".")[1][:3]
                    writer.writerow(
                        [sample_metainfo['Sample'], sample_metainfo['Description'],sample_metainfo['System'], sample_metainfo['crRNA'],sample_metainfo['Candidate alignment reads'],sample_metainfo['spike-in'],
                        sample_metainfo['mapped'], (sample_metainfo['Candidate alignment reads'] - sample_metainfo['mapped']), 
                        sample_metainfo['yield'].split(".")[0]+"."+sample_metainfo['yield'].split(".")[1][:3],
                        sample_metainfo['ontarget'], sample_metainfo['offtarget'], sample_metainfo['specificity'], sample_metainfo['unique_sites'],
                        sample_metainfo['RL'],sample_metainfo['dual_reads'],sample_metainfo['misaligned'],sample_metainfo['T7 insertions'],sample_metainfo['self_target'],sample_metainfo['mismatched_alignment'],
                        sample_metainfo['filtered_total'],sample_metainfo['filtered_on'],sample_metainfo['filtered_off']])
                    if sample_metainfo['Weblogo'] == 'Yes':
                        if sample_metainfo['unique_sites']!=0:
                            print("\n\nGenerating Fasta for Weblogo for Sample {}\n\n{} bases on each side are taken\n\n".format(sample_metainfo['Sample'],params['weblogo_flankseq']))
                            weblogo(sample_metainfo, alignment_algo, UMI_plotting=True)
                        else: print("\n\nNot generating weblogo for {}, no mapped integration sites\n\n".format(sample_metainfo['Sample']))
                    max_norm_values.append(sample_metainfo['max_bin'])
                    for item in max_norm_values:
                        if item*max_spikein>max_norm_value: max_norm_value = item*max_spikein
                    bins_dict[sample_metainfo['Sample']]=sample_metainfo['plot_bins']


print("\n----\n\nAnalyzed in {} minutes and {} seconds ".format((time.time() - start_time)//60,
    ((time.time() - start_time)%60)//1)
)
