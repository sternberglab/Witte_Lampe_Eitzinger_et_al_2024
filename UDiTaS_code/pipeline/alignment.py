import subprocess
import os
import glob

from parameters import params, index_genomes, index_folder

def generate_index_path(alignment_algo, index_folder, index_genomes, genome):
    if alignment_algo == 'bowtie2':
        index = index_folder + index_genomes[genome] + '.index'
        return [index,index]
    elif alignment_algo == 'bwa':
        index = index_folder + index_genomes[genome]
        index_check = index_folder + index_genomes[genome] + '.bwt'
        return [index_check,index]

# def align(alignment_algo, index_path, metainfo, , is_crRNA = False):
#     if not is_crRNA:
#     else: 
#         crRNA_path = 'Bowtie_input/crRNA_fasta_files/' + metainfo['Sample'] + '.fasta'

def check_index(alignment_algo, index_folder, index_genomes, genome,cores):
    if alignment_algo not in ['bwa','bowtie2']:
        print('\n\n*******\n\nNo appropriate alignment algorithm input!! please check alignemnt_algo variable in paramters.py \n\n*******\n\n')
        return False

    index_path = generate_index_path(alignment_algo, index_folder, index_genomes, genome)
    files = glob.glob(index_path[0] + '*')
    if not files:
        print("\n\n\n{} genome not indexed for {} algorithm\n\n\n".format(genome, alignment_algo))
        if alignment_algo == 'bwa':
            print(index_folder)
            os.chdir(index_folder)
            subprocess.run('bwa index {}'.format(index_genomes[genome]),shell = True)
        elif alignment_algo == 'bowtie2':
            os.chdir(index_folder)
            subprocess.run('bowtie2-build --threads {} {} {}.index'.format(cores,index_genomes[genome],index_genomes[genome]),shell = True)
        os.chdir(params['Directory'])
    else: print("\n{} genome indexed for {} algorithm".format(genome, alignment_algo))

    return True

def align(indexpath,metainfo,cores,alignment_algo,pairing):
    if not params['soft clipping'] and alignment_algo == 'bowtie2': soft_clip = " --end-to-end "
    else: soft_clip = " "
    if pairing:
        

        if params['TagTn_UMIs']:
            os.mkdir("UMI_stats/{}/".format(metainfo['Sample']))
            prededup = '{}/{}_output/{}_pre_deduplicated.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
        input_sample_1 = '{}/alignment_input/{}_input_1.fastq'.format(params['Prefix'],metainfo['Sample'])
        input_sample_2 = '{}/alignment_input/{}_input_2.fastq'.format(params['Prefix'],metainfo['Sample'])
    

    elif not pairing: input_sample = '{}/alignment_input/{}_input.fastq'.format(params['Prefix'],metainfo['Sample'])
    output_sample = '{}/{}_output/{}.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
    if len(metainfo['crRNA']) >0:
        bowtie_index = index_folder + index_genomes[metainfo['Genome']] + '.index'
        if metainfo['Genome'] == 'pace_ecoli': kmister = "-k 2"
        else: kmister = ""
        crRNA_path = params['Prefix']+'/alignment_input/crRNA_fasta_files/' + metainfo['Sample'] + '.fasta'
        crRNA_output = '{}/{}_output/'.format(params['Prefix'],alignment_algo) + metainfo['Sample'] + '_crRNA.bam'
        subprocess.run('bowtie2 -x {} -f {} -p {}{}{}--very-sensitive --quiet | samtools sort -o {} ; samtools index {}'.format(bowtie_index, crRNA_path, cores, params['bowtie2_randomness'], kmister, crRNA_output, crRNA_output), shell=True)
    if metainfo['off_targets']=='Yes':
        bowtie_index = index_folder + index_genomes[metainfo['Genome']] + '.index'
        if metainfo['Genome'] == 'pace_ecoli': bowtie_index = index_folder + index_genomes['ecoli'] + '.index'
        offtarget_path = params['Prefix']+'/alignment_input/offtarget_fasta_files/' + metainfo['Sample'] + '.fasta'
        offtarget_output = '{}/{}_output/'.format(params['Prefix'],alignment_algo) + metainfo['Sample'] + '_offtargets.bam'
        subprocess.run('bowtie2 -x {} -f {} -p {}{}--very-sensitive --quiet | samtools sort -o {} ; samtools index {}'.format(bowtie_index, offtarget_path, cores, params['bowtie2_randomness'], offtarget_output, offtarget_output), shell=True)
    if alignment_algo == 'bowtie2':
        

        if pairing and params['TagTn_UMIs']:
            subprocess.run('bowtie2 -x {} -q -1 {} -2 {} -p {}{}--very-sensitive --no-mixed --no-discordant --quiet {} | samtools sort -o {} ; samtools index {}'.format(indexpath, input_sample_1, input_sample_2,cores,soft_clip,params['bowtie2_randomness'],prededup, prededup), shell=True)
            subprocess.run('umi_tools dedup -I {} -L UMI_stats/{}/{}_log.txt --output-stats=UMI_stats/{}/{}_dedupstats -S {}'.format(prededup,metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],output_sample), shell=True)
        
        if pairing and params['pDonor_UMIs']:
            output_UMI = '{}/{}_output/{}_UMI_processed.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
            output_UMI_group = '{}/{}_output/{}_UMI_grouped.tsv'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
            output_sample_grouped = '{}/{}_output/{}_UMI_grouped.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])

            subprocess.run('bowtie2 -x {} -q -1 {} -2 {} -p {}{}--very-sensitive --no-mixed --no-discordant --quiet {} | samtools sort -o {} ; samtools index {}'.format(
                indexpath, input_sample_1, input_sample_2,cores,soft_clip,params['bowtie2_randomness'],output_sample, output_sample), shell=True)

            subprocess.run('umi_tools group -I {} --group-out={} --output-bam --paired -S{} ; samtools index {}'.format(output_sample,output_UMI_group,output_sample_grouped,output_sample_grouped,output_sample_grouped), shell=True)

            subprocess.run('umi_tools dedup -I {} --paired -S{}'.format(output_sample_grouped,output_UMI), shell=True)



            # subprocess.run('umi_tools dedup -I {} -L UMI_stats/{}/{}_log.txt --output-stats=UMI_stats/{}/{}_dedupstats -S {}'.format(
            #     output_sample,metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],ouput_UMI), shell=True)


        elif pairing: subprocess.run('bowtie2 -x {} -q -1 {} -2 {} -p {}{}{}--very-sensitive --no-mixed --no-discordant --quiet | samtools sort -o {} ; samtools index {}'.format(indexpath, input_sample_1, input_sample_2,cores,soft_clip,params['bowtie2_randomness'],output_sample, output_sample), shell=True)
        elif not pairing: subprocess.run('bowtie2 -x {} -q {} -p {}{}{}--very-sensitive --quiet | samtools sort -o {} ; samtools index {}'.format(indexpath, input_sample,cores,soft_clip,params['bowtie2_randomness'],output_sample, output_sample), shell=True)
    
    elif alignment_algo == 'bwa':
        

        if pairing and params['TagTn_UMIs']:
            output_UMI = '{}/{}_output/{}_UMI_processed.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])
            bwaalign = subprocess.run('bwa mem -t {} -v 3 {} {} {} | samtools sort -o {} ; samtools index {}'.format(cores, indexpath, input_sample_1, input_sample_2, prededup, prededup), shell=True)
            subprocess.run('umi_tools dedup -I {} -L UMI_stats/{}/{}_log.txt --output-stats=UMI_stats/{}/{}_dedupstats -S {}'.format(prededup,metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],output_sample), shell=True)
        
        if pairing and params['pDonor_UMIs']:
            output_UMI = '{}/{}_output/{}_UMI_processed.bam'.format(params['Prefix'],alignment_algo, metainfo['Sample'])

            bwaalign = subprocess.run('bwa mem -t {} -v 3 {} {} {} | samtools sort -o {} ; samtools index {}'.format(
                cores, indexpath, input_sample_1, input_sample_2, output_sample, output_sample), shell=True)
            
            subprocess.run('umi_tools group -I {} --group-out={} --output-bam --paired -S{} ; samtools index {}'.format(output_sample,output_UMI_group,output_sample_grouped,output_sample_grouped,output_sample_grouped), shell=True)

            subprocess.run('umi_tools dedup -I {} --paired -S{}'.format(output_sample_grouped,output_UMI), shell=True)

            # subprocess.run('umi_tools dedup -I {} -L UMI_stats/{}/{}_log.txt --output-stats=UMI_stats/{}/{}_dedupstats -S{}'.format(
            #     output_sample,metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],metainfo['Sample'],ouput_UMI), shell=True)


        elif pairing: bwaalign = subprocess.run('bwa mem -t {} -v 3 {} {} {} | samtools sort -o {} ; samtools index {}'.format(cores, indexpath, input_sample_1, input_sample_2, output_sample, output_sample), shell=True)
        elif not pairing: bwaalign = subprocess.run('bwa mem -t {} -v 3 {} {} | samtools sort -o {} ; samtools index {}'.format(cores, indexpath, input_sample, output_sample, output_sample), shell=True)


