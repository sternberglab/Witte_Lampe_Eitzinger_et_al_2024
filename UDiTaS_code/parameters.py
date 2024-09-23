import multiprocessing
import psutil

'''

PROCESSING ORDER FOR SAMPLE FASTQ FILES

This is currently not in use yet, but the idea is to just choose the default order based on your protein of interest. 


A GOOD COMMAND FOR CHECKING PLASMID–MAPPING READS: 


bbmap.sh in=[READS] ref=[REF] nodisk out=[OUT.bam] bamscript=bs.sh; sh bs.sh



'''

processing_order = { 
	
	'qScore' : 1,

	'Tn_end_trim' : 2,

	'pDonor_filter' : 3,

	'spikein_filter' : 4,

	'minicirlce_finder' : 5,

	'self_target' : 6,

	'plasmid_finder' : 7

}


ec2_path = ''

params = {

	### DIRECTORY MUST END WITH /

	# READ LOCATION SHOULD BE 'Reads' folder within working directory

	'Root_path': '/home/ec2-user', # this does not end with a /, and still not used

	'Code_path': '/home/ec2-user/GDL_TagTn_AWS',

	'info_file':"20240813_TagTn_UMI.csv",

	'Directory': '/home/ec2-user/xSL0324_pooled_reads/',

	'sequencing_date':'20240715',    # MUST BE YYYYMMDD format for typical sample processing to match up
	##### UPDATED TO NOT INCLUDE FOR POOLED SEQUENCES

	'Prefix':'xSL0324-pooled', # what you want the processing and alignment statistics to start with

	'zoom_value': 1.0,
	
	'is_paired_end':True,
	
	'Run_processing':True, # in case you have processing and just want to redo bam/plotting etc
	'Run_alignment':True, # in case you just want to process/trim reads for other purposes, turn this to false and the program will stop before performing bowtie alingment and plotting
	'Run_plotting':True,

	'KEEP_INTERMEDIATES':False, #  Whether intermediate FastQ and BBDuk files are kept, can help for troubleshooting, typically turned off

	'minimum_trimmed_length':20,
	'minimum_insertion_reads': 2, # number of reads at any insertion site to be counted as a "real" insertion for plotting and on/off anlaysis in the "filtered" datasets
	
	'TagTn_UMIs':False,
	'pDonor_UMIs':True,
	'UMI_length':'NNNNNNNNNN',

	'spikein_norm':False,

	'soft clipping':False, # allow soft clipping of sequencing reads (can make for sloppier alignment, typically suggest false)

	'Alignment_hamming': [True,3], # MAXIMUM amount of mismatches in the R1 that you tolerate, currently ignoring R2 cause it's far cleaner, but easily something adjustable

	#Name of info file, make sure it is IN THE WORKING DIRECTORY ABOVE

	'weblogo_flankseq':20, # MUST BE A MULTIPLE OF 5!

	'log':False, #plot on log scale

	'logscale':10, #plot on log scale

	'mismatch_tolerance':'3', # must be between 0 and 3, mismatch tolerance during sample processing

	# Typically 1 mismatch is best, but sometimes cranking to 3 on low quality runs is valuable

	'REP_strand':False,

	'plasmids':{

# T7-based BL21 plasmids: 
		# 'pEffector':{
		# 	'QCas':[1581,653,'#b4c7e3'],
		# 	'TnsABC':[6602,10365,'#f59d87']
		# },
		# 'pDonor':{
		# 	'RE':[1391,1465,'#102e59'],
		# 	'LE':[6261,6405,'#102e59']
		# },

# S2060 time-course Plasmids:
		# 'pEffector':{
		# 	'QCas':[768,5888,'#b4c7e3'],
		# 	'RE':[6034,6108,'#102e59'],
		# 	'LE':[6907,7051,'#102e59']
		# },
		# 'pTnsABC':{
		# 	'TnsABC':[1423,5186,'#f59d87']
		# },

# TnpA pEffector
		# 'pTnpA': {
		# 	'TnpA':[1071,1679,'#18735D']
		# }
	},
	'self-target_mask':[False,825,925,'pEffector'],
	'donor-end_mask':[1381,1390,'pDonor'], # PACE campain plasmid

	# 'donor-end_mask':[6025,6035,'pEffector'], # Timecourse plasmid

	'plotting_colors':{
		'bars':'#b45c4b', #CAST
		# 'bars':'#348665', #IStron
		'crRNA':'#9b2740',
		'offtargets':'#E2A28B',
		'T7':'#b6c4a2',
		'REP':'#800000',
		'REP_reverse':'#808080',
		'REP_forward':'#800000'
	}, # color to change for graphs

	'plotting_positions':{
		#'T7':365183
	}, # positions for plotting

	'font_sizes':{
	'title':7,
	'ax_label':6,
	'tick_label':6
	},
	'axis-width':.5,
	'fig_dim':2,
	'fig_ratio':3,



	'bowtie2_randomness': '',

	'chromatin': False,
	'chromatin_overlay': True,

	'bin_size': 100000,
	'sigma': 50,
	'bedgraphfile': '',
	'chromatin_color':'#348665',
	'chromatin_width':.5
}

UMI_params = {
	# How many reads any UMI needs to be counted as a real read for now, leave at 1 

	'min_UMI_read_count':2,

	# stagger length of primers, trims all reads to make sure UMI amplicons are trimmed to same length, 
	# helps ensure "UMIs" compare same amplicon lengths (otherwise not same amplicon and impacts UMI program)
	# essentially, # of cycles in a run - stagger primer length (so 142 cycle with a max stagger primer length of 5 = 137)

	# Updating code to have this trimmed down on a tiny amount of UMI loss, like 400 out of 79,000 UMIs that were discareded were then saved
	#	ofc could have been a low-depth run in general on UMIs

	'UMI_location':'id',


	## must be either "index" or "id", determines how UMI is fetched

	'max_length':130,
}

## *** if doing ecoli AND plasmid -- combine them into a single FASTA file, 'plasmid' is for if that's the only "type" of substrate used in experiment (changes plotting parameters, namely binning) ***

## *** the dictionary Key value (aka the left side, in ALL LOWERCASE), must match the genome you put in your input csv!

index_genomes = {
	'mg1655':'MG1655.fna',
	#'pace_ecoli':'tnpa_ecoli.fna',
	# 'pace_ecoli':'pace_ecoli.fna',
	# 'pace_ecoli':'pace_S2060.fna',
	# 'pace_ecoli':'pace_S2060_timecourse.fna', #LacZ region from DH10B manually removed
	#'human':'GRCh38_genome.fna',
	'ecoli':'BL21_DE3.fna',
	'human-noy':'GRCh38_genome_noy.fna', #y chromosome removed, for HEK293T cells
	'human-noy_aavs1':'GRCh38_genome_noy_AAVS1.fna', #y chromosome removed, for HEK293T cells, with HINDIII scar at AAVS1 knockin, for uSL0023 mapping, use crRNA TTATATTCCCAGGGCCGGTTAATGTGGCTCTG
	'Transfected_plasmids':'Liu_pSLs.fasta',
	'metagenome':'yeast.fna'
}

intermediate_path = "{}/Intermediates/".format(params['Prefix'])

#QScore filter (filtered on average qscore in read)

qscore = 20

# Working Directory, READS MUST BE IN HERE, will be where all folders are made

# Index path and files for genomes

# INDEX PATH MUST END WITH / 
# OTHERWISE, CODE WILL FAIL AT SPECIFIC POINTS

index_folder = "/home/ec2-user/GDL_TagTn_AWS/indexes/"

# Alignment algo ** right now, it is either 'bowtie2' or 'bwa', NO OTHER OPTIONS!!

## strongly recommended to use bowtie2

alignment_algo = 'bowtie2'

# Modify target site duplication if needed (really no reason to tinker with this)
# TSD only used for CAST, not for TnpA for now

TSD = 5

#  edge with for plotting, something to still play around with, not sure how ideal this is just yet

width = .75

# number of cores to use for bowtie, can manually change cores_to_use if needed, but otherwise don't touch
memory = psutil.virtual_memory()
total_memory = memory.total / (1024 ** 3)
available_memory = int(total_memory - 1)

cores = multiprocessing.cpu_count()
cores_to_use = max(1, cores-1)
