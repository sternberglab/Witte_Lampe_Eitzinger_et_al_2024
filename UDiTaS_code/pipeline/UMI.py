#!/usr/bin/env python3
# coding: utf-8

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import csv
import subprocess
import os
import matplotlib.pyplot as plt
import itertools

def read_process(read1,UMI_location,pairing,read2 = False,reader=([],[])):
	if not pairing:
		if UMI_location == 'id':
			UMI_seq_record = read1.id[-10:]

		elif UMI_location == 'index':
			UMI_seq_record = read1.description[-19:-9]

		return {'UMI':UMI_seq_record,
		'Sequence':str(read1.seq)}
	if pairing:

		if UMI_location == 'id':
			UMI_seq_record = read1.id[-10:]
			read1.id = read1.id[0:-11]+"_{}".format(UMI_seq_record)
			read2.id = read2.id[0:-11]+"_{}".format(UMI_seq_record)

		elif UMI_location == 'index':
			UMI_seq_record = read1.description[-19:-9]
			read1.id = read1.id + '_{}'.format(UMI_seq_record)
			read2.id = read2.id + '_{}'.format(UMI_seq_record)

		read1.description = ''
		read2.description = ''
		#print(read1)
		#print(read2)

		reader[0].append(read1)
		reader[1].append(read2)
		return {'UMI':UMI_seq_record
		# 'Sequence1':str(read1.seq),
		# 'Sequence2':str(read2.seq)
		}

def generate_UMI_plots(metainfo,UMI_counter,min_read_count,read_proportions,pairing):
	plt.rcParams['svg.fonttype'] = 'none'  # important so that text stays as intact characters in the output
	plt.rcParams['font.sans-serif'] = "Arial"
	plt.rcParams['font.family'] = "sans-serif"
	plt.rcParams['axes.spines.right']=False
	plt.rcParams['axes.spines.top']=False
	plt.rcParams['axes.spines.left']=True
	plt.rcParams['axes.spines.bottom']=True
	plt.rcParams['axes.linewidth']=.5
	binning1 = max(UMI_counter.values())//2
	if binning1 == 0:
		binning1 = 1
	if not pairing:
		fig, (ax1, ax2) = plt.subplots(2,1, dpi=600, figsize=(2,4),tight_layout=True)
		ax1.hist(UMI_counter.values(),bins=binning1, color="#9b2740",log=True)
		ax1.set_title('{}'.format(metainfo['Sample']), fontsize=6)
		ax1.text(0.5, 1.025, "{} unique UMIs".format(len(UMI_counter.keys())), fontsize=5, transform=ax1.transAxes,
		        horizontalalignment='center', verticalalignment='top')
		ax1.set_xlabel('Reads per UMI', fontsize=5)
		ax1.set_ylabel('# of UMIs', fontsize=5)
		ax1.tick_params(axis='both', which='major', labelsize=5, width=.5, length=3, direction='out',color='black',bottom=True)

		ax2.hist(read_proportions,bins=20, color="#9b2740",log=True)
		ax2.set_title('{}'.format(metainfo['Sample']), fontsize=6)
		ax2.set_xlabel('Proportion of unique sequences to total reads', fontsize=5)
		ax2.set_ylabel('# of UMIs', fontsize=5)
		ax2.tick_params(axis='both', which='major', labelsize=5, width=.5, length=3, direction='out',color='black',bottom=True)

	elif pairing:
		fig, (ax1) = plt.subplots(1,1, dpi=600, figsize=(2,2),tight_layout=True)
		ax1.hist(UMI_counter.values(),bins=binning1, color="#9b2740",log=True)
		ax1.set_title('{}'.format(metainfo['Sample']), fontsize=6)
		ax1.text(0.5, 1.025, "{} unique UMIs".format(len(UMI_counter.keys())), fontsize=5, transform=ax1.transAxes,
		        horizontalalignment='center', verticalalignment='top')
		ax1.set_xlabel('Reads per UMI', fontsize=5)
		ax1.set_ylabel('# of UMIs', fontsize=5)
		ax1.tick_params(axis='both', which='major', labelsize=5, width=.5, length=3, direction='out',color='black',bottom=True)

	plt.savefig("UMI_stats/{}_UMI.png".format(metainfo['Sample']), format='png')        
	plt.close(fig)

def UMI_output(metainfo,UMI_dictionary,UMI_counter,min_read_count,pairing,reader = ([],[])):
	output_name = metainfo['filename']
	if not pairing:
		output_file= 'Reads/' + output_name + '_UMI.fastq' # filtered file name
		output = []
		read_proportions = []
		discarded_UMI = 0
		N_discard = 0
		with open('UMI_stats/{}_UMI_statistics.csv'.format(metainfo['Sample']), 'w', newline='') as stats_file: # Open stats file for writing
		    writer = csv.writer(stats_file)
		    writer.writerow(['UMI','Reads','Unique sequences']) # Write header row
		    for UMI in list(UMI_dictionary.keys()):
		    	proportion = len(UMI_dictionary[UMI].values())/UMI_counter[UMI]
		    	read_proportions.append(proportion)
		    	writer.writerow([UMI,UMI_counter[UMI],len(UMI_dictionary[UMI])])
		    	if "N" in UMI:
		    		discarded_UMI +=1
		    		N_discard +=1
		    		continue
		    	if UMI_counter[UMI] != 1 and max(UMI_dictionary[UMI].values()) ==1:
		    		discarded_UMI += 1
		    		continue
		    	if UMI_counter[UMI] < min_read_count:
		    		discarded_UMI +=1
		    		continue
		    	sequence = max(UMI_dictionary[UMI])
		    	seq_quality = [40]*len(sequence)
		    	out_seq = SeqRecord(Seq(sequence), id = str(UMI), description = '',
		    		letter_annotations={'phred_quality': seq_quality})
		    	output.append(out_seq)
		with open(output_file, "w") as out_handle:
			SeqIO.write(output,out_handle,'fastq')
		print("\nSample {} had a total of {} UMIs and {} were discarded, {} of which were due to Ns in UMI\n\n".format(metainfo['Sample'],len(UMI_counter.keys()),discarded_UMI,N_discard))

	elif pairing:
		read_proportions = []
		output_1 = 'Reads/' + output_name + '_UMI_1.fastq'
		output_2 = 'Reads/' + output_name + '_UMI_2.fastq'

		with open(output_1, "w") as out_handle:
			SeqIO.write(reader[0],out_handle,'fastq')

		with open(output_2, "w") as out_handle:
			SeqIO.write(reader[1],out_handle,'fastq')

	generate_UMI_plots(metainfo,UMI_counter,min_read_count, read_proportions, pairing)


def UMI(metainfo,UMI_location,min_read_count,pairing):
	input_name = metainfo['filename']
	if not pairing:
		input_file='Reads/'+input_name+'.fastq'
		input_data = SeqIO.parse(input_file,'fastq')
		UMI_dictionary = {}
		UMI_counter = {}
		for input_read in input_data:
			read_data = read_process(input_read,UMI_location,pairing)
			if read_data['UMI'] in UMI_dictionary:
				UMI_counter[read_data['UMI']] +=1
				if read_data['Sequence'] in UMI_dictionary[read_data['UMI']]:
					UMI_dictionary[read_data['UMI']][read_data['Sequence']] +=1
				else:
					UMI_dictionary[read_data['UMI']][read_data['Sequence']] =1
			else:
				UMI_counter[read_data['UMI']] = 1
				UMI_dictionary[read_data['UMI']] = {}
				UMI_dictionary[read_data['UMI']][read_data['Sequence']] =1
		UMI_output(metainfo,UMI_dictionary,UMI_counter,min_read_count,pairing)
	elif pairing:
		input_file_1='Reads/'+input_name+'_1.fastq'
		input_file_2='Reads/'+input_name+'_2.fastq'
		UMI_dictionary = {}
		UMI_counter = {}
		input_data_1 = SeqIO.parse(input_file_1,'fastq')
		input_data_2 = SeqIO.parse(input_file_2,'fastq')
		reader=([],[])
		for (input_r1, input_r2) in zip(input_data_1,input_data_2):
			read_data = read_process(input_r1, UMI_location,pairing,input_r2,reader)
			if read_data['UMI'] in UMI_counter:
				UMI_counter[read_data['UMI']] +=1
			else:
				UMI_counter[read_data['UMI']] = 1

		UMI_output(metainfo,UMI_dictionary,UMI_counter,min_read_count,pairing,reader)

