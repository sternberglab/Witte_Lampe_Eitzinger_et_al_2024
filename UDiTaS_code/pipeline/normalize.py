import pysam
import matplotlib.pyplot as plt
import os
import numpy as np
import csv
from Bio import SeqIO
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from parameters import params, index_folder, index_genomes
from pipeline.alignment import generate_index_path
from pipeline.plotting import check_read,get_orientation_distance,plot_features,crRNA_map,crRNA_map_nonhuman,off_target_map,setup_yaxis,plotbins,process_insertion,analyze_insertion

def setup_plot(axs,maxy,metainfo,plot_crRNA,offset,x_marks,x_labels,data):
    axs.xaxis.set_ticks_position('none')
    axs.set_ylim(-.12*maxy,maxy)
    axs.set_xlim(0-(5*10**7),offset+(5*10**7))
    axs.set_ylabel("Percent of {}s".format(data), fontsize=5)
    axs.tick_params(axis='both', which='major', labelsize=5, width=.75, length=3, direction='out',color='black',bottom=True)
    axs.tick_params(axis='x', rotation=90)
    axs.set_xticks(x_marks)
    axs.set_xticklabels(x_labels)
    axs.set_xlabel("Chromosome", fontsize=5)
    axs.set_title(metainfo['Description'], fontsize=6)
    axs.margins(1,1)
    axs.set_yticks(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.set_yticklabels(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.spines['left'].set_bounds(0, maxy)
    axs.spines['left'].set_position(('data',(0-5*10**7)))
    axs.spines['bottom'].set_position('zero')
    axs.tick_params(axis='x', which='minor', width=.5, length=2, direction='out',color='black',bottom=True)
    if metainfo['System'] == 'CAST' and plot_crRNA[0] == True:
        crRNA_bin = 500 + ((plot_crRNA[1] // 1000) * 1000)
        axs.scatter(plot_crRNA[1], -.05 * maxy, marker="^", c="#9b2740", s=45, edgecolor='black', linewidth = .25)

def setup_contig_plot(axs,maxy,metainfo,plot_crRNA,offset,data,PLASMID=False): # fix to have 4.5 vs 4 handled
    offset_scientific = str("{0:.2E}".format(float(offset)))
    exp = offset_scientific[offset_scientific.index("+")+1:]
    if exp[0]=="0": exp = exp[1:]
    base = offset_scientific[:offset_scientific.index(".")]
    base_tens = float(offset_scientific[:4])
    axs.xaxis.set_ticks_position('none')
    axs.spines['left'].set_position(('data',(0-5*offset*10**(-3))))
    axs.spines['bottom'].set_position('zero')
    axs.set_ylabel("Percent of {}s".format(data), fontsize=5)
    axs.set_ylim(-.12*maxy,maxy)
    axs.tick_params(axis='both', which='major', labelsize=5, width=.75, length=3, direction='out',color='black',bottom=True)
    axs.set_xlabel("Genome Position ($10^{}$ bp)".format(exp), fontsize=5)
    if PLASMID: axs.set_xlabel("Plasmid Position ($10^{}$ bp)".format(exp), fontsize=5)
    axs.set_title(metainfo['Description'], fontsize=6)
    axs.set_xlim((0-5*offset*10**(-3)),offset+5*offset*10**(-3))
    axs.set_xticks(np.arange(0,base_tens*10**int(exp),.5*10**int(exp)))
    axs.set_xticklabels(np.arange(0,base_tens,.5))
    axs.margins(1,1)
    axs.xaxis.set_minor_locator(MultipleLocator(.25*10**int(exp)))
    axs.tick_params(axis='x', which='minor', width=.5, length=2, direction='out',color='black',bottom=True)
    axs.set_title(metainfo['Description'], fontsize=6)
    axs.margins(1,1)
    axs.set_yticks(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.set_yticklabels(np.append(np.arange(0,maxy,maxy/2),maxy))
    axs.spines['left'].set_bounds(0, maxy)
    axs.spines['left'].set_position(('data',(0-5*offset*10**(-3))))
    axs.spines['bottom'].set_position('zero')
    if metainfo['System'] == 'CAST' and plot_crRNA[0] == True:
        if metainfo['Genome'] in ['pace_ecoli','ecoli']:
            for crRNA_hit in plot_crRNA[1]:
                if crRNA_hit[0] == 'BL21': 
                    crRNA_bin = 500 + ((crRNA_hit[1] // 1000) * 1000)
                    axs.scatter(crRNA_hit[1], -.05 * maxy, marker="^", c="#9b2740", s=45, edgecolor='black', linewidth = .25)
                else: continue
        else:
            crRNA_bin = 500 + ((plot_crRNA[1] // 1000) * 1000)
            axs.scatter(plot_crRNA[1], -.05 * maxy, marker="^", c="#9b2740", s=45, edgecolor='black', linewidth = .25)

def percentage(plot_bins_percentages,plot_bins,metainfo,max_spikein):
    for item in list(plot_bins.keys()):
        plot_bins_percentages[item] = {}
        for location in plot_bins[item]: plot_bins_percentages[item][location]=max_spikein*plot_bins[item][location]/int(metainfo['spike-in'])

def plotnormalizefig(metainfo, plot_bins,width,data,max_spikein,max_norm_value,alignment_algo):
    maxy = max_norm_value
    label_range = [0,maxy/2,maxy]
    labels = [0,50,100]
    plot_bins_percentages = {}
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['axes.spines.right']=False
    plt.rcParams['axes.spines.top']=False
    plt.rcParams['axes.spines.left']=True
    plt.rcParams['axes.spines.bottom']=True
    plt.rcParams['axes.linewidth']=.75

    bamfile = pysam.AlignmentFile("{}/{}_output/{}.bam".format(params['Prefix'],alignment_algo,metainfo['Sample']), "rb")

    ref_names = bamfile.references
    ref_lengths = bamfile.lengths

    if metainfo['Genome'] in ['pace_ecoli','ecoli']:
        plot_crRNA = crRNA_map_nonhuman(metainfo,alignment_algo)
        percentage(plot_bins_percentages,plot_bins,metainfo,max_spikein)
        max_chromosomes = []
        if 'BL21' in list(plot_bins_percentages.keys()): max_chromosomes.append(max(plot_bins['BL21'].values()))
        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4), dpi=1000, tight_layout=True)
        for plot in [[ax1,maxy],[ax2,.005*maxy]]: setup_contig_plot(plot[0],plot[1],metainfo,plot_crRNA,4.6*10**6,data)
        ax1.scatter(751500, -.05*maxy, marker="^",c='#b6c4a2',s=45, edgecolor='black',linewidth=.5)
        ax1.set_yticks(label_range)
        ax1.set_yticklabels(labels)
        ax2.set_yticks([0,params['zoom_value']*.01*maxy])
        ax2.set_yticklabels([0,params['zoom_value']*.01])
        if metainfo['System'] == 'CAST' and plot_crRNA[0] == True:
            for crRNA_hit in plot_crRNA[1]:
                if crRNA_hit[0] == 'BL21':
                    crRNA_bin = 500 + ((crRNA_hit[1] // 1000) * 1000)
                    for ax in [ax1,ax2]: ax.scatter(crRNA_hit[1], -.05 * maxy, marker="^", c="#9b2740", s=45, edgecolor='black', linewidth = .55)
        ax2.scatter(751500, -.05*params['zoom_value']*.01*maxy, marker="^",c='#b6c4a2',s=45, edgecolor='black',linewidth=.5)
        if 'BL21' in list(plot_bins_percentages.keys()): plotbins('BL21',plot_bins_percentages,ax1,ax2,width)
        for ax in [ax1,ax2]: ax.set_ylabel("Normalized reads", fontsize=5)
        plt.savefig("{}/Output_figures/Normalized/{}_normalized.pdf".format(params['Prefix'],metainfo['Sample'],metainfo['Sample']), format='pdf')
        plt.close(fig)
        for plasmid in list(plot_bins.keys()):
            if plasmid == "BL21": continue
            for ref_name, ref_length in zip(ref_names, ref_lengths):
                if plasmid == ref_name: x_range = ref_length
            fig, ax = plt.subplots(1, figsize=(8, 2), dpi=1000, tight_layout=True)
            setup_contig_plot(ax,maxy,metainfo,plot_crRNA,x_range,data,PLASMID=True)
            ax.set_title("{} – {}".format(metainfo['Description'],plasmid), fontsize=6)
            if len(list(params['plasmids'].keys()))>0 and plasmid in list(params['plasmids'].keys()): plot_features(ax,params['plasmids'],plasmid)
            setup_yaxis(0,maxy,ax)
            ax.set_yticks(label_range)
            ax.set_yticklabels(labels)
            for location in plot_bins_percentages[plasmid].keys():
                if location>params['self-target_mask'][1] and location<params['self-target_mask'][2] and plasmid == params['self-target_mask'][3]: continue
                elif location>=params['donor-end_mask'][0] and location<=params['donor-end_mask'][1] and plasmid == params['donor-end_mask'][2]: continue
                else: ax.bar(location,plot_bins_percentages[plasmid][location],edgecolor='#b45c4b',linewidth=width)
            if plot_crRNA[0]:
                for crRNA_hit in plot_crRNA[1]:
                    if crRNA_hit[0] == plasmid: ax.scatter(crRNA_hit[4], -.05*maxy, marker="^",c='#9b2740',s=45, edgecolor='black',linewidth=.5)
            ax.set_ylabel("Normalized reads", fontsize=5)
            plt.savefig("{}/Output_figures/Normalized/{}_{}_normalized.pdf".format(params['Prefix'],metainfo['Sample'],plasmid), format='pdf')
            plt.close(fig)

    elif metainfo['Genome'] in ['human','human-noy','metagenome']:
        chromosomes = ['NC_000001.11','NC_000002.12','NC_000003.12','NC_000004.12','NC_000005.10','NC_000006.12','NC_000007.14','NC_000008.11','NC_000009.12','NC_000010.11','NC_000011.10','NC_000012.12','NC_000013.11','NC_000014.9','NC_000015.10','NC_000016.10','NC_000017.11','NC_000018.10','NC_000019.10','NC_000020.11','NC_000021.9','NC_000022.11','NC_000023.11','NC_000024.10']
        chromosome_labels = ['Chromosome 1','Chromosome 2','Chromosome 3','Chromosome 4','Chromosome 5','Chromosome 6','Chromosome 7','Chromosome 8','Chromosome 9','Chromosome 10','Chromosome 11','Chromosome 12','Chromosome 13','Chromosome 14','Chromosome 15','Chromosome 16','Chromosome 17','Chromosome 18','Chromosome 19','Chromosome 20','Chromosome 21','Chromosome 22','Chromosome X','Chromosome Y']
        if metainfo['Genome'] == 'human-noy':
            chromosomes = chromosomes[:-1]
            chromosome_labels = chromosome_labels[:-1]
        x_marks = []
        offset_lengths = {}
        offset = 0
        for ref_name, ref_length in zip(ref_names, ref_lengths):
            offset_lengths[ref_name] = offset
            if ref_name in chromosomes: x_marks.append(offset)
            offset += ref_length
        plot_crRNA = crRNA_map(metainfo,offset_lengths,alignment_algo)
        offtargets = off_target_map(metainfo,offset_lengths,alignment_algo)
        percentage(plot_bins_percentages,plot_bins,metainfo,max_spikein)
        x_labels=[]
        for i in chromosome_labels: x_labels.append(i.split(" ")[1])
        fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4), dpi=1000, tight_layout=True)
        max_chromosomes = []
        for item in plot_bins_percentages.keys(): max_chromosomes.append(max(plot_bins_percentages[item].values()))
        if metainfo['Genome'] in ['human','human-noy']: 
            for plot in [[ax1,maxy],[ax2,params['zoom_value']*.01*maxy]]: setup_plot(plot[0],plot[1],metainfo,plot_crRNA,offset,x_marks,x_labels,data)
        else: 
            for plot in [[ax1,maxy],[ax2,params['zoom_value']*.01*maxy]]:setup_contig_plot(plot[0],plot[1],metainfo,plot_crRNA,offset,data)
        if metainfo['System'] == 'CAST' and metainfo['off_targets']=='Yes' and offtargets:
            for item in offtargets.keys():
                ax1.scatter(offtargets[item][0], -.05 * maxy, marker="^", c="#E2A28B", s=30, edgecolor='black', linewidth = .25)
                ax2.scatter(offtargets[item][0], -.05 * maxy*params['zoom_value']*.01, marker="^", c="#E2A28B", s=30, edgecolor='black', linewidth = .25)
        ax1.set_yticks(label_range)
        ax1.set_yticklabels(labels)
        ax2.set_yticks([0,params['zoom_value']*.01*maxy])
        ax2.set_yticklabels([0,params['zoom_value']])
        for key in plot_bins_percentages.keys(): plotbins(key,plot_bins_percentages,ax1,ax2,width)
        for ax in [ax1,ax2]: ax.set_ylabel("Normalized reads", fontsize=5)
        plt.savefig("{}/Output_figures/Normalized/{}_normalized.pdf".format(params['Prefix'],metainfo['Sample'],metainfo['Sample']), format='pdf')
        plt.close(fig)

