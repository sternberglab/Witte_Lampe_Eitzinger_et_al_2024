import pysam
import matplotlib.pyplot as plt
import os
import matplotlib.font_manager
import numpy as np
import csv
from Bio import SeqIO
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from parameters import params, index_folder, index_genomes, UMI_params
from pipeline.alignment import generate_index_path
from pipeline.chromatin import *

def filtered_specificity(metainfo,crRNA_info,counter,read_positions,filtered_read_positions,filtered_read_positions_percentages,filtered_plot_bins,filtered_plot_bins_percentages,data):
    for key in read_positions.keys():
        if data == 'UMI': min_insertions = 1
        else: min_insertions = params['minimum_insertion_reads']
        filtered_read_positions[key] = {k:v for k,v in read_positions[key].items() if v >= min_insertions}
        if len(filtered_read_positions[key])==0:
            del filtered_read_positions[key]
    for key in filtered_read_positions.keys():
        counter['filtered_total']+= sum(filtered_read_positions[key].values())
        filtered_read_positions_percentages[key]={}
        filtered_plot_bins[key]={}
        filtered_plot_bins_percentages[key]={}
    for key in filtered_read_positions.keys():
        for start_pos in filtered_read_positions[key].keys():
            filtered_read_positions_percentages[key][start_pos]=0
            counter['filtered_unique_sites']+=1
            if crRNA_info[2] and crRNA_info[1] < start_pos and start_pos-crRNA_info[1]<=100: counter['filtered_on'] += filtered_read_positions[key][start_pos] 
            elif not crRNA_info[2] and crRNA_info[1] > start_pos and crRNA_info[1]-start_pos<=100: counter['filtered_on'] +=filtered_read_positions[key][start_pos] 
            else: counter['filtered_off'] += filtered_read_positions[key][start_pos]
            binned_site = 500 + ((start_pos // 1000) * 1000)
            if binned_site not in filtered_plot_bins[key].keys():
                filtered_plot_bins[key][binned_site] = filtered_read_positions[key][start_pos]
                filtered_plot_bins_percentages[key][binned_site] = 0
            else: filtered_plot_bins[key][binned_site] += filtered_read_positions[key][start_pos] 
    percentage(filtered_read_positions_percentages,filtered_read_positions,filtered_plot_bins_percentages,filtered_plot_bins,counter['filtered_total'])

def check_read(read,pairing,unmapped_fasta,mapped_fasta,counter,alignment_algo,start_pos,data,umi_groups,UMI_plotting):
    mismatch = False
    UMI_read_depth = True
    if read.is_read2:
        if not read.is_unmapped:
            for item in read.get_tags():
                if item[0] == 'NM': counter['mismatch2_dist'].append(item[1])
        return False
    elif pairing and read.template_length > 1200:
        counter['misaligned'] += 1
        return False
    elif read.is_unmapped:
        unmapped_fasta.append(SeqRecord(Seq(read.get_forward_sequence()),id=read.query_name,description = ""))
        return False
    elif read.is_secondary:
        counter['secondary']+=1
        return False
    tags = read.get_tags()
    if alignment_algo == 'bowtie2' and tags[1][0] == 'XS' and tags[0][1] == tags[1][1]: counter['dual_reads'] += 1
    for item in tags:
        if item[0] == 'NM':
            counter['mismatch_dist'].append(item[1])
            if int(item[1])>params['Alignment_hamming'][1] and params['Alignment_hamming'][0]:
                counter['mismatched_alignment']+=1
                mismatch = True
        if item[0] == 'UG' and UMI_plotting:
            if umi_groups[item[1]] >= UMI_params['min_UMI_read_count']: UMI_read_depth = True
            elif umi_groups[item[1]] < UMI_params['min_UMI_read_count']: UMI_read_depth =  False
    if not UMI_read_depth: return False
    if mismatch: return False
    if alignment_algo == 'bwa' and len(tags) > 4:
        if tags[4][0] == 'XA': counter['dual_reads'] += 1
    counter['mapped'] += 1
    mapped_fasta.append(SeqRecord(Seq(read.get_forward_sequence()),id = read.query_name,description=read.reference_name+"_"+str(start_pos)))
    return True

def get_orientation_distance(crRNA_info,start_pos,read_strand,TSD):
    if crRNA_info[2] and crRNA_info[1] < start_pos:
        if read_strand:
            distance_from_crRNA = start_pos - crRNA_info[1] - crRNA_info[3] + TSD
            orientation = 'LR'
        else:
            distance_from_crRNA = start_pos - crRNA_info[1] - crRNA_info[3]
            orientation = 'RL'   
    elif not crRNA_info[2] and crRNA_info[1] > start_pos:
        if read_strand:
            distance_from_crRNA = crRNA_info[1]-start_pos
            orientation = 'RL' 
        else:
            distance_from_crRNA = crRNA_info[1] - start_pos + TSD
            orientation = 'LR'
    else:
        distance_from_crRNA = 1000
        orientation = 'N/A'
    return [distance_from_crRNA, orientation]

def update_insertion_profile(insertion_info,counter,RL_sample_distribution,LR_sample_distribution,ontarget_fasta,offtarget_fasta,read):
    if insertion_info[0]<=100:
        ontarget_fasta.append(SeqRecord(Seq(read.get_forward_sequence()),id=read.query_name,description = ""))
        counter['ontarget'] += 1
        if insertion_info[1] == 'RL':
            if insertion_info[0] not in RL_sample_distribution: RL_sample_distribution[insertion_info[0]] = 1
            else: RL_sample_distribution[insertion_info[0]] += 1
            counter['RL_counter'] +=1
        elif insertion_info[1] == 'LR':
            if insertion_info[0] not in LR_sample_distribution: LR_sample_distribution[insertion_info[0]] = 1
            else: LR_sample_distribution[insertion_info[0]] += 1
            counter['LR_counter'] +=1
    else:
        offtarget_fasta.append(SeqRecord(Seq(read.get_forward_sequence()),id=read.query_name,description = ""))
        counter['offtarget'] += 1

def analyze_insertion(plot_crRNA,metainfo,start_pos,TSD,counter,RL_sample_distribution,LR_sample_distribution,ontarget_fasta,offtarget_fasta,read):
    if plot_crRNA[0]:
        if metainfo['Genome'] in ['pace_ecoli','ecoli']:
            for crRNA_hit in plot_crRNA[1]:
                if crRNA_hit[0]==read.reference_name:
                    insertion_info = get_orientation_distance(crRNA_hit,start_pos,read.is_forward,TSD)
                    break
                else: insertion_info = [1000,'N/A']
        else: insertion_info = get_orientation_distance(plot_crRNA,start_pos,read.is_forward,TSD)
        update_insertion_profile(insertion_info,counter,RL_sample_distribution,LR_sample_distribution,ontarget_fasta,offtarget_fasta,read)
    elif not plot_crRNA[0]:
        counter['ontarget'] = False
        counter['offtarget'] = False
        counter['RL_counter'] = False
        counter['LR_counter'] = False

def plot_features(ax,plasmids,plasmid):
    for sub in plasmids[str(plasmid)].keys(): ax.axvspan(plasmids[str(plasmid)][sub][0], plasmids[str(plasmid)][sub][1], ymin=.107, facecolor=plasmids[str(plasmid)][sub][2], edgecolor='None',alpha=0.25,label=sub)

def crRNA_map(metainfo,offset_lengths, alignment_algo):
    if metainfo['System'] != 'CAST': 
        print("\n\nNot a CAST sample, no crRNA\n\n")
        return [False,False,False,False] 
    try:
        crRNA_bamfile = pysam.AlignmentFile("{}/{}_output/{}_crRNA.bam".format(params['Prefix'],alignment_algo, metainfo['Sample']), "rb")
        for read in crRNA_bamfile:
            # Check if the read is aligned
            if read.is_unmapped:
                print("\nNo crRNA alignment occured! \n")
                return [False, False, False]
                continue
            crRNA_location = offset_lengths[read.reference_name] + read.reference_start
            crRNA_strand = read.is_forward
            crRNA_length = read.query_length
            return [True, crRNA_location, crRNA_strand, crRNA_length]
    except Exception as e:
        print('\nNo crRNA FASTA file for {}! \n'.format(metainfo['Sample']))
        return [False,False,False,False]  

def crRNA_map_nonhuman(metainfo, alignment_algo):
    if metainfo['System'] != 'CAST': 
        print("\nNot a CAST sample, no crRNA\n")
        return [False,False,False,False] 
    try:
        crRNA_bamfile = pysam.AlignmentFile("{}/{}_output/{}_crRNA.bam".format(params['Prefix'],alignment_algo, metainfo['Sample']), "rb")
        crRNA_hits = []
        for read in crRNA_bamfile:
            # Check if the read is aligned
            if read.is_unmapped:
                print("\nThe crRNA failed to align! \n")
                return [False, False, False]
                continue
            crRNA_contig = read.reference_name
            crRNA_location = read.reference_start
            crRNA_strand = read.is_forward
            crRNA_length = read.query_length
            if crRNA_strand: crRNA_plotting = crRNA_location
            else: crRNA_plotting = crRNA_location + crRNA_length
            crRNA_hits.append([crRNA_contig,crRNA_location, crRNA_strand, crRNA_length, crRNA_plotting])
        return [True, crRNA_hits]
    except Exception as e:
        print('\nNo crRNA FASTA file for {}! \n'.format(metainfo['Sample']))
        return [False,False,False,False]  

def off_target_map(metainfo,offset_lengths, alignment_algo):
    try:
        offtarget_alignments = {}
        counter = 1
        offtarget_bamfile = pysam.AlignmentFile('{}/{}_output/'.format(params['Prefix'],alignment_algo) + metainfo['Sample'] + '_offtargets.bam', "rb")
        for read in offtarget_bamfile:
            # Check if the read is aligned
            if read.is_unmapped:
                print("\nNo off-target alignment occured! \n")
                return [False, False, False]
                continue
            offtarget_location = offset_lengths[read.reference_name] + read.reference_start
            offtarget_alignments[counter] = [offtarget_location, read.is_forward, read.query_length]
            counter += 1
        return offtarget_alignments
    except Exception as e:
        print('\nNo off-target FASTA file for {}! \n'.format(metainfo['Sample']))
        return False 

def setup_plot(axs,maxy,metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter,chromatin=False):
    axs.xaxis.set_ticks_position('none')
    if not chromatin: 
        axs.set_ylim(-.12*maxy,maxy)
        axs.spines['left'].set_bounds(0, maxy)
    if chromatin: 
        if not params['chromatin_overlay']:
            axs.set_ylim(1,maxy)
            axs.spines['left'].set_bounds(1, maxy)
        else: 
            axs.set_ylim(.9925,maxy*1.025)
    axs.set_xlim(0-(5*10**7),offset+(5*10**7))
    axs.set_ylabel("Percent of {}s".format(data), fontsize=params['font_sizes']['ax_label'])
    axs.tick_params(axis='both', which='major', labelsize=params['font_sizes']['tick_label'], width=params['axis-width'], length=3, direction='out',color='black',bottom=True)
    axs.tick_params(axis='x', rotation=90)
    axs.set_xticks(x_marks)
    axs.set_xticklabels(x_labels)
    axs.set_xlabel("Chromosome", fontsize=params['font_sizes']['ax_label'])
    if not chromatin: axs.set_title(metainfo['Description']+", {} Mapped {}s".format(counter['mapped'],data), fontsize=params['font_sizes']['title'])
    axs.margins(1,1)
    if not chromatin: 
        axs.set_yticks(np.append(np.arange(0,maxy,maxy/2),maxy))
        axs.set_yticklabels(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.spines['left'].set_position(('data',(0-5*10**7)))
    if not chromatin: axs.spines['bottom'].set_position('zero')
    axs.tick_params(axis='x', which='minor', width=params['axis-width']*.5, length=2, direction='out',color='black',bottom=True)
    if metainfo['System'] == 'CAST' and plot_crRNA[0] == True:
        crRNA_bin = 500 + ((plot_crRNA[1] // 1000) * 1000)
        axs.scatter(plot_crRNA[1], -.05 * maxy, marker="^", c=params['plotting_colors']['crRNA'], s=45, edgecolor='black', linewidth = .25)


def setup_contig_plot(axs,maxy,metainfo,plot_crRNA,offset,data,counter,PLASMID=False): # fix to have 4.5 vs 4 handled
    offset_scientific = str("{0:.2E}".format(float(offset)))
    exp = offset_scientific[offset_scientific.index("+")+1:]
    if exp[0]=="0": exp = exp[1:]
    base = offset_scientific[:offset_scientific.index(".")]
    base_tens = float(offset_scientific[:4])
    axs.xaxis.set_ticks_position('none')
    axs.spines['left'].set_position(('data',(0-5*offset*10**(-3))))
    axs.spines['bottom'].set_position('zero')
    axs.set_ylabel("Percent of {}s".format(data), fontsize=params['font_sizes']['ax_label'])
    axs.set_ylim(-.12*maxy,maxy)
    axs.tick_params(axis='both', which='major', labelsize=params['font_sizes']['tick_label'], width=params['axis-width'], length=3, direction='out',color='black',bottom=True)
    axs.set_xlabel("Genome Position ($10^{}$ bp)".format(exp), fontsize=params['font_sizes']['ax_label'])
    if PLASMID: axs.set_xlabel("Plasmid Position ($10^{}$ bp)".format(exp), fontsize=params['font_sizes']['ax_label'])
    axs.set_title(metainfo['Description'], fontsize=params['font_sizes']['title'])
    axs.set_xlim((0-5*offset*10**(-3)),offset+5*offset*10**(-3))
    axs.set_xticks(np.arange(0,base_tens*10**int(exp),.5*10**int(exp)))
    axs.set_xticklabels(np.arange(0,base_tens,.5))
    axs.margins(1,1)
    axs.xaxis.set_minor_locator(MultipleLocator(.25*10**int(exp)))
    axs.tick_params(axis='x', which='minor', width=params['axis-width']*.5, length=2, direction='out',color='black',bottom=True)
    axs.set_title(metainfo['Description']+", {} Mapped {}s".format(counter['mapped'],data), fontsize=params['font_sizes']['title'])
    axs.margins(1,1)
    axs.set_yticks(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.set_yticklabels(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.spines['left'].set_bounds(0, maxy)
    if metainfo['System'] == 'CAST' and plot_crRNA[0]:
        if metainfo['Genome'] in ['pace_ecoli','ecoli']:
            for crRNA_hit in plot_crRNA[1]:
                if crRNA_hit[0] in ['BL21','DH10B','MG1655']: 
                    crRNA_bin = 500 + ((crRNA_hit[1] // 1000) * 1000)
                    axs.scatter(crRNA_hit[1], -.05 * maxy, marker="^", c=params['plotting_colors']['crRNA'], s=45, edgecolor='black', linewidth = .25)
                else: continue
        else:
            crRNA_bin = 500 + ((plot_crRNA[1] // 1000) * 1000)
            axs.scatter(plot_crRNA[1], -.05 * maxy, marker="^", c=params['plotting_colors']['crRNA'], s=45, edgecolor='black', linewidth = .25)

def plot_distribution(RL_sample_distribution,LR_sample_distribution,metainfo,data,output_path):
    sample_values = []
    sample_distribution = []
    for item in list(RL_sample_distribution.keys()):
        if item < 35 or item > 65: continue
        else:
            sample_distribution.append(item)
            sample_values.append(RL_sample_distribution[item])
    for item in list(LR_sample_distribution.keys()):
        if item < 35 or item > 65: continue
        else:
            sample_distribution.append(item)
            sample_values.append(LR_sample_distribution[item])
    if len(sample_distribution) != 0:
        maxx = ((max(sample_distribution)+9) // 5)* 5
        minx = (min(sample_distribution) // 5)* 5
        maxy = (max(sample_values))
    else:
        maxx = 60
        minx = 40
        maxy = 100
    fig = plt.figure(figsize=(2, 2), dpi=600, tight_layout=True)
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_title('{}'.format(metainfo['Description']),fontsize=params['font_sizes']['title'])
    ax1.tick_params(axis='both', which='major', labelsize=params['font_sizes']['tick_label'], width=params['axis-width'], length=3, direction='out',color='black',bottom=True)
    ax1.set_xlabel('Distance (bp)', fontsize=params['font_sizes']['ax_label'])
    ax1.set_ylabel('{}s'.format(data), fontsize=params['font_sizes']['ax_label'])
    ax1.set_xlim(minx,maxx)
    setup_yaxis(0,maxy,ax1)
    ax1.set_ylim(0,maxy)
    ax1.margins(1,1)
    ax1.set_xticks(np.append(np.arange(minx,maxx,(maxx-minx)/(2*(maxx-minx)/10)),maxx))
    if len(RL_sample_distribution.keys()) != 0: plt.bar(RL_sample_distribution.keys(), RL_sample_distribution.values(), color='#b45c4b', label='T-RL')
    else: plt.bar(0,0, color='#b45c4b', label='T-RL')
    if len(LR_sample_distribution.keys()) != 0: plt.bar(LR_sample_distribution.keys(), LR_sample_distribution.values(), facecolor='None', edgecolor = 'black', linewidth=.5, label='T-LR')
    else: plt.bar(0,0, facecolor='None', edgecolor = 'black', linewidth=.5, label='T-LR')            
    handles, labels = ax1.get_legend_handles_labels()
    ax1_legend = ax1.twinx()
    ax1_legend.spines[['top','right','left','bottom']].set_visible(False)
    ax1_legend.axis('off')
    legend = ax1_legend.legend(handles, labels, fontsize = params['font_sizes']['ax_label'], handlelength=1,frameon=False, loc='center left', bbox_to_anchor=(1.05, .9))
    legend.get_frame().set_linewidth(0)
    plt.savefig("{}/{}/Distributions/{}_Distribution.pdf".format(params['Prefix'],output_path,metainfo['Sample']), format='pdf')        
    plt.close(fig)

def accuracy_output(plot_crRNA, metainfo, counter, alignment_algo, ontarget_fasta, offtarget_fasta,unmapped_fasta,mapped_fasta):
    SeqIO.write(unmapped_fasta, '{}/{}_output/alignment_outputs/{}_unmapped.fasta'.format(params['Prefix'],alignment_algo, metainfo['Sample']), "fasta")
    SeqIO.write(mapped_fasta, '{}/{}_output/alignment_outputs/{}_mapped.fasta'.format(params['Prefix'],alignment_algo, metainfo['Sample']), "fasta")
    if not plot_crRNA[0]:
        counter['ontarget'] = False
        counter['offtarget'] = False
        counter['RL_counter'] = False
        counter['LR_counter'] = False
    if plot_crRNA[0]:
        print("\n------\n\n{} total ontarget reads for Sample {}\n".format(counter['ontarget'], metainfo['Sample']))
        print("\n{} total offtarget reads for Sample {}\n\n------\n".format(counter['offtarget'], metainfo['Sample']))
        SeqIO.write(ontarget_fasta, '{}/{}_output/alignment_outputs/{}_ontarget.fasta'.format(params['Prefix'],alignment_algo, metainfo['Sample']), "fasta")
        SeqIO.write(offtarget_fasta, '{}/{}_output/alignment_outputs/{}_offtarget.fasta'.format(params['Prefix'],alignment_algo, metainfo['Sample']), "fasta")
    print("\n------\n\n{} total reads for Sample {}\n".format(metainfo['Candidate alignment reads'], metainfo['Sample']))
    print("\n{} total mapped reads for Sample {}\n\n------\n".format(counter['mapped'], metainfo['Sample']))

def process_insertion(metainfo,read,read_positions,csv_positions,read_orientations,plot_bins,read_positions_percentages,plot_bins_percentages,counter,start_pos,csv_pos):
    if read.reference_name not in read_positions.keys():
        csv_positions[read.reference_name]={}
        read_orientations[read.reference_name]={}
        plot_bins[read.reference_name]={}
        plot_bins_percentages[read.reference_name]={}
        read_positions_percentages[read.reference_name]={}
        read_positions[read.reference_name]={}
    if start_pos not in read_positions[read.reference_name].keys():
        read_positions_percentages[read.reference_name][start_pos]=0
        csv_positions[read.reference_name][csv_pos] = 1
        read_orientations[read.reference_name][csv_pos] = [0,0]
        read_positions[read.reference_name][start_pos] = 1
        counter['unique_sites'] +=1
        if read.reference_name in ['BL21','DH10B','MG1655'] or metainfo['Genome'] in ['human','human-noy',"human-noy_aavs1"]: binned_site = 500 + ((start_pos // 1000) * 1000)
        else: binned_site = 5 + ((start_pos // 10) * 10)
        if binned_site not in plot_bins[read.reference_name].keys():
            plot_bins[read.reference_name][binned_site] = 1
            plot_bins_percentages[read.reference_name][binned_site] = 0
        else: plot_bins[read.reference_name][binned_site] += 1
    else:
        read_positions[read.reference_name][start_pos] += 1
        csv_positions[read.reference_name][csv_pos] += 1
        if read.reference_name in ['DH10B','BL21','MG1655'] or metainfo['Genome'] in ['human','human-noy','human-noy_aavs1']: binned_site = 500 + ((start_pos // 1000) * 1000)
        else: binned_site = 5 + ((start_pos // 10) * 10)
        plot_bins[read.reference_name][binned_site]+=1
    if read.is_forward: read_orientations[read.reference_name][csv_pos][1]+=1
    else: read_orientations[read.reference_name][csv_pos][0]+=1

def setup_yaxis(miny,maxy,axi):
    axi.set_ylim(-.12*maxy,maxy)
    axi.set_yticks([miny,maxy//2,maxy])
    axi.set_yticklabels([miny,maxy//2,maxy])
    axi.spines['left'].set_bounds(miny, maxy)

def plotbins(key,plot_bins_percentages,ax1,ax2,width):
    for location in plot_bins_percentages[key]:
        ax1.bar(location,plot_bins_percentages[key][location], edgecolor=params['plotting_colors']['bars'],linewidth=width, zorder=2)
        ax2.bar(location,plot_bins_percentages[key][location], edgecolor=params['plotting_colors']['bars'],linewidth=width, zorder=2)

def percentage(read_positions_percentages,read_positions,plot_bins_percentages,plot_bins,counter):
    for item in read_positions_percentages.keys():
        for location in read_positions_percentages[item]: read_positions_percentages[item][location]=100*read_positions[item][location]/counter
    for item in plot_bins_percentages.keys():
        for location in plot_bins_percentages[item]: plot_bins_percentages[item][location]=100*plot_bins[item][location]/counter

def get_max_norm(counter,plot_bins,metainfo):
    temp_max = 0
    for section in plot_bins:
        for region in list(plot_bins[section].keys()):
            if plot_bins[section][region]/int(metainfo['spike-in'])>temp_max: temp_max = plot_bins[section][region]/int(metainfo['spike-in'])
    counter['max_bin'] = temp_max

def plot_chromatin(params,axs,max_atac,metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter,positions,smoothed_coverage):
    setup_plot(axs,max_atac,metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter,chromatin=True)
    axs.plot(positions, smoothed_coverage, color=params['chromatin_color'],linewidth=params['chromatin_width'],zorder=1,alpha=.6)
    axs.tick_params(axis='y', left=False, right=False)
    axs.set_yticklabels([])
    axs.set_ylabel("")
    axs.spines['left'].set_visible(False)  # Hide the y-axis spine
    if params['chromatin_overlay']: axs.spines['bottom'].set_visible(False)


def plotfig(metainfo, TSD, alignment_algo,width,data,pairing,UMI_plotting=False):

    if UMI_plotting: output_path = "UMI_output_figures"
    else: output_path = "Output_figures"
    umi_groups = {}
    if UMI_plotting: 
        UMI_tsv = '{}/{}_output/{}_UMI_grouped.tsv'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
        with open(UMI_tsv, 'r') as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')
            for row in reader:
                umi_group = int(row['unique_id'])
                read_count = int(row['final_umi_count'])
                umi_groups[umi_group] = read_count

    counter = {
    'misaligned':0,
    'mapped' : 0,
    'ontarget': 0,
    'offtarget': 0,
    'RL_counter' : 0,
    'LR_counter' : 0,
    'unique_sites': 0,
    'dual_reads':0,
    'T7 insertions':0,
    'self_target':0,
    'secondary':0,
    'mismatch_dist':[],
    'mismatch2_dist':[],
    'mismatched_alignment':0,
    'max_bin':0,
    'resid_donor':0,
    'plot_bins':{},
    'filtered_on':0,
    'filtered_off':0,
    'filtered_total':0,
    'filtered_unique_sites':0
    }
    RL_sample_distribution = {}
    LR_sample_distribution = {}
    read_positions = {}
    read_positions_percentages = {}
    filtered_read_positions_percentages = {}
    plot_bins = {}
    plot_bins_percentages = {}
    filtered_read_positions = {}
    filtered_plot_bins = {}
    filtered_plot_bins_percentages = {}
    read_orientations = {}
    csv_positions = {}
    unmapped_fasta = []
    mapped_fasta = []
    ontarget_fasta = []
    offtarget_fasta = []
    # plt.rcParams['font'] = "sans-serif"
    plt.rcParams['font.family'] = "DejaVu Sans"
    plt.rcParams['axes.spines.right']=False
    plt.rcParams['axes.spines.top']=False
    plt.rcParams['axes.spines.left']=True
    plt.rcParams['axes.spines.bottom']=True
    plt.rcParams['axes.linewidth']=params['axis-width']

    bamfile = pysam.AlignmentFile("{}/{}_output/{}.bam".format(params['Prefix'],alignment_algo,metainfo['Sample']), "rb")
    if UMI_plotting: bamfile = pysam.AlignmentFile("{}/{}_output/{}_UMI_processed.bam".format(params['Prefix'],alignment_algo,metainfo['Sample']), "rb")

    ref_names = bamfile.references
    ref_lengths = bamfile.lengths

    if metainfo['Genome'] in ['human','human-noy','metagenome','human-noy_aavs1']:
        chromosomes = ['NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10','NC_000006.12','NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10']
        labels = ['Chromosome 1','Chromosome 2','Chromosome 3','Chromosome 4','Chromosome 5','Chromosome 6','Chromosome 7','Chromosome 8','Chromosome 9','Chromosome 10','Chromosome 11','Chromosome 12','Chromosome 13','Chromosome 14','Chromosome 15','Chromosome 16','Chromosome 17','Chromosome 18','Chromosome 19','Chromosome 20','Chromosome 21','Chromosome 22','Chromosome X','Chromosome Y']
        if metainfo['Genome'] == 'human-noy' or metainfo['Genome'] == 'human-noy_aavs1':
            chromosomes = chromosomes[:-1]
            labels = labels[:-1]
        x_marks = []
        offset_lengths = {}
        offset = 0
        for ref_name, ref_length in zip(ref_names, ref_lengths):
            offset_lengths[ref_name] = offset
            if ref_name in chromosomes: x_marks.append(offset)
            offset += ref_length
        plot_crRNA = crRNA_map(metainfo,offset_lengths,alignment_algo)
        offtargets = off_target_map(metainfo,offset_lengths,alignment_algo)
        for read in bamfile:
            # Check if the read is aligned
            if read.is_unmapped: start_pos = 0
            elif read.is_forward: start_pos = offset_lengths[read.reference_name] + read.reference_start
            else: start_pos = offset_lengths[read.reference_name] + read.reference_start + read.query_length
            checks = check_read(read,pairing,unmapped_fasta,mapped_fasta,counter,alignment_algo,start_pos,data,umi_groups,UMI_plotting)
            if not checks: continue
            csv_pos = start_pos - offset_lengths[read.reference_name]
            analyze_insertion(plot_crRNA,metainfo,start_pos,TSD,counter,RL_sample_distribution,LR_sample_distribution,ontarget_fasta,offtarget_fasta,read)
            process_insertion(metainfo,read,read_positions,csv_positions,read_orientations,plot_bins,read_positions_percentages,plot_bins_percentages,counter,start_pos,csv_pos)
        accuracy_output(plot_crRNA, metainfo, counter, alignment_algo, ontarget_fasta, offtarget_fasta,unmapped_fasta,mapped_fasta)
        percentage(read_positions_percentages,read_positions,plot_bins_percentages,plot_bins,counter['mapped'])

        filtered_specificity(metainfo,plot_crRNA,counter,read_positions,filtered_read_positions,filtered_read_positions_percentages,filtered_plot_bins,filtered_plot_bins_percentages,data)

        x_labels=[]

        if params['chromatin']: 
            positions, smoothed_coverage = generate_smoothed_coverage(params['bedgraphfile'], bin_size=params['bin_size'], sigma=params['sigma'])

        for i in labels: x_labels.append(i.split(" ")[1])
        fig, (ax1, ax2) = plt.subplots(2, figsize=(params['fig_dim']*params['fig_ratio'],params['fig_dim']*.5*params['fig_ratio']), dpi=1000, tight_layout=True)
        if metainfo['Genome'] in ['human','human-noy','human-noy_aavs1']: 
            for plot in [[ax1,100],[ax2,.5]]: setup_plot(plot[0],plot[1],metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter)
        else: 
            for plot in [[ax1,100],[ax2,.5]]:setup_contig_plot(plot[0],plot[1],metainfo,plot_crRNA,offset,data,counter)
        if metainfo['System'] == 'CAST' and metainfo['off_targets']=='Yes' and offtargets:
            for item in offtargets.keys():
                ax1.scatter(offtargets[item][0], -.05 * 100, marker="^", c=params['plotting_colors']['offtargets'], s=30, edgecolor='black', linewidth = .25)
                ax2.scatter(offtargets[item][0], -.05 * .5, marker="^", c=params['plotting_colors']['offtargets'], s=30, edgecolor='black', linewidth = .25)
        for key in plot_bins_percentages.keys(): plotbins(key,plot_bins_percentages,ax1,ax2,width)
        if plot_crRNA[0] and (counter['ontarget']+counter['offtarget']) != 0:
            specificity = str(100*counter['ontarget']/(counter['ontarget']+counter['offtarget']))
            specificity = specificity.split(".")[0]+"."+specificity.split(".")[1][:3]
            ax1.text(0.5, 1.025, specificity + "%" + " on-target {}s".format(data), fontsize=params['font_sizes']['ax_label'], transform=ax1.transAxes, horizontalalignment='center', verticalalignment='top')
        plt.savefig("{}/{}/Percent/{}_percent.pdf".format(params['Prefix'],output_path,metainfo['Sample'],metainfo['Sample']), format='pdf')
        plt.close(fig)
        if plot_crRNA[0] and counter['ontarget']>=1: plot_distribution(RL_sample_distribution,LR_sample_distribution,metainfo,data,output_path)


        if params['chromatin']:
            if params['chromatin_overlay']: 
                fig2, (ax3) = plt.subplots(1, figsize=(params['fig_dim']*params['fig_ratio'],params['fig_dim']), dpi=1000, tight_layout=True)
                ax6 = ax3.twinx()
            else: fig2, (ax3,ax6) = plt.subplots(2, figsize=(params['fig_dim']*params['fig_ratio'],1.25*params['fig_dim']), gridspec_kw={'height_ratios': [4, 1]}, dpi=1000, tight_layout=True)        
        else: fig2, (ax3) = plt.subplots(1, figsize=(params['fig_dim']*params['fig_ratio'],params['fig_dim']), dpi=1000, tight_layout=True)
        
        max_chromosomes = []
        for item in read_positions.keys(): max_chromosomes.append(max(plot_bins[item].values()))
        if len(max_chromosomes)==0: maxy = 1
        else: maxy = max(max_chromosomes)
        if metainfo['Genome'] in ['human','human-noy','human-noy_aavs1']: setup_plot(ax3,maxy,metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter)
        else: setup_contig_plot(ax3,maxy,metainfo,plot_crRNA,offset,data,counter)
        if metainfo['System'] == 'CAST' and metainfo['off_targets']=='Yes':
            if offtargets:
                for item in offtargets.keys(): ax3.scatter(offtargets[item][0], -.05 * maxy, marker="^", c=params['plotting_colors']['offtargets'], s=30, edgecolor='black', linewidth = .25)
        for key in plot_bins.keys():
            for location in plot_bins[key]: ax3.bar(location,plot_bins[key][location], edgecolor=params['plotting_colors']['bars'],linewidth=width, zorder=2)
        ax3.set_ylabel("Number of {}s".format(data), fontsize=params['font_sizes']['ax_label'])
        if plot_crRNA[0] and (counter['ontarget']+counter['offtarget']) != 0:
            specificity = str(100*counter['ontarget']/(counter['ontarget']+counter['offtarget']))
            specificity = specificity.split(".")[0]+"."+specificity.split(".")[1][:3]
            ax3.text(0.5, 1.025, specificity + "%" + " on-target {}s".format(data), fontsize=params['font_sizes']['ax_label'], transform=ax3.transAxes,
                    horizontalalignment='center', verticalalignment='top')
        if params['chromatin']: 
            plot_chromatin(params,ax6,get_max_y_ATAC_value(smoothed_coverage),metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter,positions,smoothed_coverage)
            if not params['chromatin_overlay']: ax6.set_xlabel("")  # Hide x-axis label


        plt.savefig("{}/{}/Raw/{}_{}s.pdf".format(params['Prefix'],output_path,metainfo['Sample'],data), format='pdf')
        plt.close(fig2)

        if params['chromatin']:
            if params['chromatin_overlay']: 
                fig3, (ax4) = plt.subplots(1, figsize=(params['fig_dim']*params['fig_ratio'],params['fig_dim']), dpi=1000, tight_layout=True)
                ax5 = ax4.twinx()
            else: fig3, (ax4,ax5) = plt.subplots(2, figsize=(params['fig_dim']*params['fig_ratio'],1.25*params['fig_dim']), gridspec_kw={'height_ratios': [4, 1]}, dpi=1000, tight_layout=True)
        
        else: fig3, (ax4) = plt.subplots(1, figsize=(params['fig_dim']*params['fig_ratio'],params['fig_dim']), dpi=1000, tight_layout=True)
        max_chromosomes = []
        for item in filtered_read_positions.keys(): max_chromosomes.append(max(filtered_plot_bins[item].values()))
        if len(max_chromosomes)==0: maxy = 1
        else: maxy = max(max_chromosomes)
        if metainfo['Genome'] in ['human','human-noy','human-noy_aavs1']: setup_plot(ax4,maxy,metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter)
        else: setup_contig_plot(ax4,maxy,metainfo,plot_crRNA,offset,data,counter)
        if metainfo['System'] == 'CAST' and metainfo['off_targets']=='Yes':
            if offtargets:
                for item in offtargets.keys(): ax4.scatter(offtargets[item][0], -.05 * maxy, marker="^", c=params['plotting_colors']['offtargets'], s=30, edgecolor='black', linewidth = .25)
        for key in filtered_plot_bins.keys():
            for location in filtered_plot_bins[key]: ax4.bar(location,filtered_plot_bins[key][location], edgecolor=params['plotting_colors']['bars'],linewidth=width, zorder=2)
        ax4.set_ylabel("Number of {}s".format(data), fontsize=params['font_sizes']['ax_label'])
        if plot_crRNA[0] and (counter['filtered_total']) != 0:
            specificity = str(100*counter['filtered_on']/counter['filtered_total'])
            specificity = specificity.split(".")[0]+"."+specificity.split(".")[1][:3]
            ax4.text(0.5, 1.025, specificity + "%" + " on-target {}s".format(data), fontsize=params['font_sizes']['ax_label'], transform=ax4.transAxes,
                    horizontalalignment='center', verticalalignment='top')
        if params['chromatin']: 
            plot_chromatin(params,ax5,get_max_y_ATAC_value(smoothed_coverage),metainfo,plot_crRNA,offset,x_marks,x_labels,data,counter,positions,smoothed_coverage)
            if not params['chromatin_overlay']: ax4.set_xlabel("")  # Hide x-axis label

        ax4.set_title(metainfo['Description']+", {} Mapped filtered {}s".format(counter['filtered_total'],data), fontsize=params['font_sizes']['title'])
        plt.savefig("{}/{}/Filtered_Raw/{}_{}s.pdf".format(params['Prefix'],output_path,metainfo['Sample'],data), format='pdf')
        plt.close(fig3)

        with open('{}/{}/CSV_Files/{}_plotting_statistics.csv'.format(params['Prefix'],output_path,metainfo['Sample']), 'w', newline='') as stats_file: # Open stats file for writing
            writer = csv.writer(stats_file)
            writer.writerow(['Site','Insertion Count','RL','LR'])
            for key in csv_positions.keys():
                for location in csv_positions[key]: writer.writerow(["{}_{}".format(key,location),csv_positions[key][location],read_orientations[key][location][0],read_orientations[key][location][1]])

        with open('{}/{}/CSV_Files/{}_alignment_statistics.csv'.format(params['Prefix'],output_path,metainfo['Sample']), 'w', newline='') as stats_file: # Open stats file for writing
            writer = csv.writer(stats_file)
            writer.writerow(["R1","R2"])
            for (item1,item2) in zip(counter['mismatch_dist'],counter['mismatch2_dist']): writer.writerow([item1,item2])
        if params['spikein_norm']: get_max_norm(counter,plot_bins,metainfo)
    counter['plot_bins'] = plot_bins
    counter['filtered_sites']=filtered_read_positions
    return counter

def gen_label(parameter):
    logo_n = parameter['weblogo_flankseq']
    logo_length = logo_n*2
    logo_label = [',']*logo_length
    iterator = 1
    for i in range(len(logo_label)):
        if i<logo_n:
            if i == (logo_n-1): logo_label[i] = "-1,"
            if i%5==0: logo_label[i] = str(int(-1*logo_n)+int(5*i/5))+","
        elif i==logo_n: logo_label[i] = "1,"
        elif i>logo_n:
            if (i+1)%5==0:
                logo_label[i] = str(int(5*iterator))+","
                if i == logo_length-1: logo_label[i] = str(int(5*iterator))
                iterator +=1
    return "".join(logo_label)

def weblogo(metainfo, alignment_algo,UMI_plotting=False):
    if UMI_plotting: output_path = "UMI_output_figures"
    else: output_path = "Output_figures"
    weblogo_label = gen_label(params)
    index_path = index_folder + index_genomes[metainfo['Genome']]
    genome = SeqIO.parse(index_path,'fasta')
    bamfile = pysam.AlignmentFile('{}/{}_output/{}.bam'.format(params['Prefix'],alignment_algo,metainfo['Sample']),'rb')
    
    ref_names = bamfile.references
    ref_lengths = bamfile.lengths

    output_handle = '{}/{}/WeblogoFastas/{}.fasta'.format(params['Prefix'],output_path,metainfo['Sample'])
    filtered_output_handle = '{}/{}/WeblogoFastas/{}_filtered.fasta'.format(params['Prefix'],output_path,metainfo['Sample'])
    genome_ids = []
    genome_chromosomes = []
    reader = []
    filtered_reader = []
    integration_sites = []
    offset_lengths = {}
    offset = 0
    for ref_name, ref_length in zip(ref_names, ref_lengths):
        offset_lengths[ref_name] = offset
        offset += ref_length
    for record in genome:
        genome_ids.append(record.name)
        genome_chromosomes.append(record)
    for read in bamfile:
        filtered_read=False
        if read.is_read2 or read.is_unmapped: continue
        if read.is_forward: start_pos = offset_lengths[read.reference_name] + read.reference_start
        else: start_pos = offset_lengths[read.reference_name] + read.reference_start + read.query_length
        if start_pos in integration_sites: continue
        integration_sites.append(start_pos)
        for key in metainfo['filtered_sites']:
            if start_pos in metainfo['filtered_sites'][key].keys():
                filtered_read = True
        if read.reference_name in genome_ids:
            for record in genome_chromosomes:
                if read.reference_name == record.name:
                    read_name = str(read.query_name)
                    if read.is_forward:
                        read_seq = record.seq[(read.reference_start-params['weblogo_flankseq']):(read.reference_start+(params['weblogo_flankseq']))]
                    else:
                        read_seq = record.seq[(read.reference_start+read.query_length-params['weblogo_flankseq']):(read.reference_start+read.query_length+params['weblogo_flankseq'])].reverse_complement()
                    if len(read_seq) < 2*params['weblogo_flankseq']: continue
                    sequence = SeqRecord(read_seq, id = read_name, description = '')
                    reader.append(sequence)
                    if filtered_read: filtered_reader.append(sequence)

    SeqIO.write(reader, output_handle, 'fasta')
    SeqIO.write(filtered_reader, filtered_output_handle, 'fasta')

    output_path_weblogo = '{}/{}/Weblogos/{}'.format(params['Prefix'], output_path,metainfo['Sample'])
    filtered_output_path = '{}/{}/Weblogos_Filtered/{}'.format(params['Prefix'],output_path,metainfo['Sample'])


    general_params = '-F pdf -c classic --errorbars NO --resolution 600 --fineprint ""' 

    subprocess.run('weblogo -f {} -o {}_revcomp.pdf {} --revcomp --annotate {} --title "{}  {} Unique insertions"'.format(output_handle, output_path_weblogo, general_params, weblogo_label,metainfo['Description'],len(reader)), shell=True)
    subprocess.run('weblogo -f {} -o {}.pdf {} --annotate {} --title "{}  {} Unique insertions"'.format(output_handle, output_path_weblogo, general_params, weblogo_label,metainfo['Description'],len(reader)), shell=True)

    subprocess.run('weblogo -f {} -o {}_revcomp.pdf {} --revcomp --annotate {} --title "{}  {} Filtered insertions"'.format(filtered_output_handle, filtered_output_path, general_params, weblogo_label,metainfo['Description'],len(filtered_reader)), shell=True)
    subprocess.run('weblogo -f {} -o {}.pdf {} --annotate {} --title "{}  {} Filtered insertions"'.format(filtered_output_handle, filtered_output_path, general_params, weblogo_label,metainfo['Description'],len(filtered_reader)), shell=True)

    if metainfo['System'] == 'CAST':
        range_logo_start = params['weblogo_flankseq']-9
        range_logo_end = params['weblogo_flankseq']+5
        commas_left = "".join([","]*(params['weblogo_flankseq']-11))
        commas_right = "".join([","]*(params['weblogo_flankseq']-6))
        TSD_label = ",-5,-4,-3,-2,-1,1,2,3,4,5,+1,+2,+3,+4,+5,"
        labels = commas_left+TSD_label+commas_right

        subprocess.run('weblogo -f {} -o {}_lowbit2.pdf {} -l {} -u {} --revcomp -S .2 --stack-width 15 --annotate {} --title "{}  {} Unique insertions .2 Bit"'.format(output_handle, output_path_weblogo, general_params, range_logo_start,range_logo_end,labels,metainfo['Description'],len(reader)), shell=True)
        subprocess.run('weblogo -f {} -o {}_lowbit5.pdf {} -l {} -u {} --revcomp -S .5 --stack-width 15 --annotate {} --title "{}  {} Unique insertions .5 Bit"'.format(output_handle, output_path_weblogo, general_params, range_logo_start,range_logo_end,labels,metainfo['Description'],len(reader)), shell=True)

        subprocess.run('weblogo -f {} -o {}_lowbit2.pdf {} -l {} -u {} --revcomp -S .2 --stack-width 15 --annotate {} --title "{}  {} Unique insertions .2 Bit"'.format(filtered_output_handle, filtered_output_path, general_params, range_logo_start,range_logo_end,labels,metainfo['Description'],len(filtered_reader)), shell=True)
        subprocess.run('weblogo -f {} -o {}_lowbit5.pdf {} -l {} -u {} --revcomp -S .5 --stack-width 15 --annotate {} --title "{}  {} Unique insertions .5 Bit"'.format(filtered_output_handle, filtered_output_path, general_params, range_logo_start,range_logo_end,labels,metainfo['Description'],len(filtered_reader)), shell=True)



