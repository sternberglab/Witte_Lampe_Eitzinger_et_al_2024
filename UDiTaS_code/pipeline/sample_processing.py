import subprocess
import os

from parameters import params, UMI_params, qscore, index_folder,index_genomes,intermediate_path,available_memory

def bbduk_func_paired(input_file,output_unmatched,filter_seq,output_matched=False):
    if len(filter_seq)<20: kmer = len(filter_seq)
    else: kmer = 20
    if not output_matched:   
        return subprocess.run(['bbduk.sh', 'in1={}_1.fastq'.format(input_file), 'in2={}_2.fastq'.format(input_file), 
            'out1={}_1.fastq out2={}_2.fastq'.format(output_unmatched,output_unmatched),  
            'k={}'.format(kmer), 'hdist='+params['mismatch_tolerance'], 'literal=' + filter_seq, 'rcomp=f', 'minlen={}'.format(params['minimum_trimmed_length']), '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
    else:
        return subprocess.run(['bbduk.sh', 'in1={}_1.fastq'.format(input_file), 'in2={}_2.fastq'.format(input_file), 
            'out1={}_1.fastq out2={}_2.fastq'.format(output_unmatched,output_unmatched), 
            'outm1={}_1.fastq outm2={}_2.fastq'.format(output_matched,output_matched), 
            'k={}'.format(kmer), 'hdist='+params['mismatch_tolerance'], 'literal=' + filter_seq, 'rcomp=f', 'minlen={}'.format(params['minimum_trimmed_length']), '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)


def bbduk_func_single(input_file,output_matched,output_unmatched,filter_seq,stats):
    print("blah")

def read_len(input_file):
    with open(input_file, 'r') as output:
        lines = output.readlines()
        return(int(len(lines) / 4))

def generate_sample_metainfo(row, sequencing_date):
    sample_metainfo = {
    "Sample": row[0],
    "Description": row[1],
    "System": row[2],
    "Genome": row[3].lower(),
    'Primer': row[4],
    "tn_end": row[5],
    "tn_end_RC": row[6],
    "crRNA": row[7],
    "pDonor_seq": row[8],
    "spikein_seq": row[9],
    "CRISPRarray": row[10],
    "Weblogo": row[11],
    "Minicircle": row[12],
    "off_targets": row[13].split(","),
    'Minicircle_counts': 'N/A',
    "filename": row[0],
    # "filename": sequencing_date + "_" + row[0],
    "Genome_mapping_reads": True
    }
    if sample_metainfo['off_targets']==['']:sample_metainfo['off_targets']=[]
    return sample_metainfo

def get_yield(inputs,output):
    yields = 100*output/inputs
    yields = str(f'{yields:.5f}')
    return yields

def umi_tools(metainfo,pairing=params['is_paired_end'],UMI_sequence=params['UMI_length']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform UMI processing for!\n")
        return [0,0,0]
    name = metainfo['Sample']
    primer = metainfo['Primer']
    UMI_seq = UMI_sequence
    if not pairing:
        
        input_file=intermediate_path + '{}/{}_filtered.fastq'.format(name,name)
        output_unmatched_file = intermediate_path + '{}/{}_primer_unmatched.fastq'.format(name,name) # Output file name for reads that don't have transposon end
        primer_file = intermediate_path + '{}/{}_primer.fastq'.format(name,name) # filtered file name
        trim_file = intermediate_path + '{}/{}_primer_trimmed.fastq'.format(name,name) # trimmed file name
        UMI_file = intermediate_path + '{}/{}_UMI_processed.fastq'.format(name,name) # trimmed file name

        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched_file, 'outm=' + primer_file,
            'k={}'.format(len(primer)), 'hdist='+params['mismatch_tolerance'], 'literal=' + primer, 'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        bbduk_process_1 = subprocess.run(['bbduk.sh', 'in=' + primer_file, 'out=' + trim_file,
            'ktrim=l', 'k={}'.format(len(primer)), 
            'literal=' + primer, 'rcomp=f', 'hdist='+params['mismatch_tolerance'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
        
        bbduk_process_2 = subprocess.run(['umi_tools', 'extract', '-I{}'.format(trim_file), '--bc-pattern={}'.format(UMI_seq),
            '--stdout={}_1.fastq'.format(UMI_file)], stderr=subprocess.PIPE, text=True)# UMI PROCESS HERE

        input_reads = read_len(input_file)
        retained_reads = read_len(UMI_file)
    
    elif pairing:

        input_file = intermediate_path + '{}/{}_filtered'.format(name,name)
        output_unmatched = intermediate_path + '{}/{}_primer_unmatched'.format(name,name)
        primer_file = intermediate_path + '{}/{}_primer'.format(name,name)
        trim_file = intermediate_path + '{}/{}_primer_trimmed'.format(name,name)
        UMI_file = intermediate_path + '{}/{}_UMI_processed'.format(name,name)

        bbduk_process = bbduk_func_paired(input_file,output_unmatched,primer,primer_file)

        bbduk_process_1 = subprocess.run(['bbduk.sh', 'in1={}_1.fastq in2={}_2.fastq'.format(primer_file,primer_file),
            'out1={}_1.fastq out2={}_2.fastq'.format(trim_file,trim_file),
            'ktrim=l', 'k={}'.format(len(primer)),
            'literal=' + primer, 'rcomp=f', 'hdist='+params['mismatch_tolerance'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        bbduk_process_2 = subprocess.run(['umi_tools', 'extract', '-I{}_1.fastq'.format(trim_file), '--bc-pattern={}'.format(UMI_seq),
            '--read2-in={}_2.fastq'.format(trim_file), '--stdout={}_1.fastq'.format(UMI_file), '--read2-out={}_2.fastq'.format(UMI_file)], stderr=subprocess.PIPE, text=True)

        input_reads = read_len('{}_1.fastq'.format(input_file))
        retained_reads = read_len('{}_1.fastq'.format(UMI_file))

    yields = get_yield(input_reads,retained_reads)

    with open(intermediate_path + '{}/{}_primer_finder_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process.stderr)
    with open(intermediate_path + '{}/{}_primer_trimmer_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process_1.stderr)
    with open(intermediate_path + '{}/{}_UMI_process_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process_2.stderr)

    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [input_reads,retained_reads, "0.000"]

    return [input_reads,retained_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3]]


def quality_filter(metainfo,pairing=params['is_paired_end']):
    input_name = metainfo['filename']
    output_name = metainfo['Sample']

    if not pairing:
        input_file='Reads/'+input_name+'.fastq'
        filtered_file = intermediate_path + '{}/{}_filtered.fastq'.format(output_name,output_name)
        file_exists = os.path.isfile(input_file)
        if file_exists: bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + filtered_file, 'maq={}'.format(qscore),'-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        else:
            print("\nNo file found for {}\n".format(input_name))
            metainfo["Genome_mapping_reads"] = False
            return [0,0,0]

        input_reads = read_len(input_file)
        retained_reads = read_len(filtered_file)
    
    elif pairing:
        input_file = 'Reads/'+input_name
        filtered_file = intermediate_path + '{}/{}_filtered'.format(output_name,output_name)

        file1_exists = os.path.isfile(input_file+'_R1.fastq')

        file2_exists = os.path.isfile(input_file+"_R2.fastq")
        
        if file1_exists and file2_exists:
            bbduk_process = subprocess.run(['bbduk.sh', 'in1={}_R1.fastq'.format(input_file), 'in2={}_R2.fastq'.format(input_file), 
                'out1={}_1.fastq out2={}_2.fastq'.format(filtered_file,filtered_file), 
                'maq={}'.format(qscore), '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
        else:
            print("\nNo paired end files found for {}\n".format(input_name))
            metainfo["Genome_mapping_reads"] = False
            return [0,0,0]

        input_reads = read_len('{}_R1.fastq'.format(input_file))
        retained_reads = read_len('{}_1.fastq'.format(filtered_file))

    yields = get_yield(input_reads,retained_reads)

    with open(intermediate_path + '{}/{}_QscoreFiltering_bbduk.txt'.format(output_name,output_name),'w') as f: f.write(bbduk_process.stderr)
    
    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [input_reads,retained_reads, "0.000"]

    return [input_reads,retained_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3]]

'''
Stagger trimming needs to be adjusted for pairing!
'''
def stagger_trim(metainfo, stagger = UMI_params['max_length']):
    input_name = metainfo['filename']
    output_name = metainfo['filename']
    primer = metainfo['Primer']
    stagger_length = stagger - len(primer)
    input_file='Reads/staggered_raw/'+input_name+'.fastq'
    staggered_file = 'Reads/' + input_name + '.fastq' # filtered file name
    bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + staggered_file, 'k={}'.format(len(primer)), 'hdist='+params['mismatch_tolerance'], 'literal=' + primer, 'ktrim=l', 'rcomp=f', 'trimq={}'.format(qscore), '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

    input_2 = 'Reads/staggered_trimmed/' + input_name + '.fastq'
    staggered_filter = 'Reads/' + output_name + '.fastq' # filtered file name
    bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_2, 'out=' + staggered_filter, 'ftr={}'.format(stagger_length), '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

def Tn_finder_trimmer(metainfo,min_length=params['minimum_trimmed_length'],pairing=params['is_paired_end'],performed_UMI=params['pDonor_UMIs']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform Tn trimming for!\n")
        return [0,0,0]
    name = metainfo['Sample']
    tn_endseq = metainfo['tn_end']
    if not pairing:
        
        input_file=intermediate_path + '{}/{}_filtered.fastq'.format(name,name)
        if performed_UMI: input_file=intermediate_path + '{}/{}_UMI_processed.fastq'.format(name,name)
        output_unmatched_file = intermediate_path + '{}/{}_unmatched.fastq'.format(name,name) # Output file name for reads that don't have transposon end
        tn_file = intermediate_path + '{}/{}_Tn.fastq'.format(name,name) # filtered file name
        trim_file = intermediate_path + '{}/{}_TnTrimmed.fastq'.format(name,name) # trimmed file name

        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched_file, 'outm=' + tn_file,
            'k={}'.format(len(tn_endseq)), 'hdist='+params['mismatch_tolerance'], 'literal=' + tn_endseq, 'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        bbduk_process_1 = subprocess.run(['bbduk.sh', 'in=' + tn_file, 'out=' + trim_file,
            'ktrim=l', 'k={}'.format(len(tn_endseq)), 
            'literal=' + tn_endseq, 'rcomp=f', 'minlen={}'.format(min_length), 'hdist='+params['mismatch_tolerance'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
        
        input_reads = read_len(input_file)
        retained_reads = read_len(trim_file)
    
    elif pairing:

        input_file = intermediate_path + '{}/{}_filtered'.format(name,name)
        if performed_UMI: input_file=intermediate_path + '{}/{}_UMI_processed'.format(name,name)
        output_unmatched = intermediate_path + '{}/{}_unmatched'.format(name,name)
        tn_file = intermediate_path + '{}/{}_Tn'.format(name,name)
        trim_file_preRC = intermediate_path + '{}/{}_TnTrimmed_noRC'.format(name,name)
        trim_file = intermediate_path + '{}/{}_TnTrimmed'.format(name,name)
        tn_rev = metainfo['tn_end_RC']

        bbduk_process = bbduk_func_paired(input_file,output_unmatched,tn_endseq,tn_file)

        bbduk_process_1 = subprocess.run(['bbduk.sh', 'in1={}_1.fastq in2={}_2.fastq'.format(tn_file,tn_file),
            'out1={}_1.fastq out2={}_2.fastq'.format(trim_file_preRC,trim_file_preRC),
            'ktrim=l', 'k={}'.format(len(tn_endseq)),
            'literal=' + tn_endseq, 'rcomp=f', 'minlen={}'.format(min_length), 'hdist='+params['mismatch_tolerance'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        bbduk_process_2 = subprocess.run(['bbduk.sh', 'in1={}_1.fastq in2={}_2.fastq'.format(trim_file_preRC,trim_file_preRC),
            'out1={}_1.fastq out2={}_2.fastq'.format(trim_file,trim_file),
            'ktrim=r', 'k={}'.format(len(tn_rev)),
            'literal=' + tn_rev, 'rcomp=f', 'minlen=10', 'hdist='+params['mismatch_tolerance'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        input_reads = read_len('{}_1.fastq'.format(input_file))
        retained_reads = read_len('{}_1.fastq'.format(trim_file))

    yields = get_yield(input_reads,retained_reads)

    with open(intermediate_path + '{}/{}_Tnfinder_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process.stderr)
    with open(intermediate_path + '{}/{}_Tnfilter_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process_1.stderr)
    if pairing:
        with open(intermediate_path + '{}/{}_TnfilterRC_bbduk.txt'.format(name,name),'w') as f: f.write(bbduk_process_2.stderr)

    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [input_reads,retained_reads, "0.000"]

    return [input_reads,retained_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3]]

def pDonor_finder(metainfo, stats,pairing=params['is_paired_end']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform Donor trimming for!\n")
        return ["N/A", "N/A", 0]
    
    input_name = metainfo['Sample']
    donor_seq = metainfo['pDonor_seq']

    if not pairing:
        input_file=intermediate_path + '{}/{}_TnTrimmed.fastq'.format(input_name,input_name)
        if metainfo['System'] == "CAST" or metainfo['System'] == "TnpA": output_unmatched = intermediate_path + '{}/{}_nodonor.fastq'.format(input_name,input_name) 
        else: output_unmatched = '{}/alignment_input/{}_input.fastq'.format(params['Prefix'],input_name) # Output file name for reads that don't have transposon end
        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched,
            'literal=' + donor_seq, 'k=15', 'hdist='+params['mismatch_tolerance'],
            'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
        retained_reads = read_len(output_unmatched)

    elif pairing:
        input_file=intermediate_path + '{}/{}_TnTrimmed'.format(input_name,input_name)
        if metainfo['System'] == "CAST" or metainfo['System'] == "TnpA": output_unmatched = intermediate_path + '{}/{}_nodonor'.format(input_name,input_name)
        else: output_unmatched = '{}/alignment_input/{}_input'.format(params['Prefix'],input_name)

        bbduk_process = bbduk_func_paired(input_file,output_unmatched,donor_seq)

        retained_reads = read_len("{}_1.fastq".format(output_unmatched))

    donor_reads = stats[1] - retained_reads

    yields = get_yield(stats[1],donor_reads)

    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [donor_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3],retained_reads]

    with open(intermediate_path + '{}/{}_donortrimmed_bbduk.txt'.format(input_name,input_name),'w') as f: f.write(bbduk_process.stderr)

    return [donor_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3],retained_reads]

def spikein_finder(metainfo,stats,pairing=params['is_paired_end']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform Spikein trimming for!\n")
        return ['NA', 'NA',0]
    input_name = metainfo['Sample']
    spikein_seq = metainfo['spikein_seq']
    if len(spikein_seq)<20: kmer = len(spikein_seq)
    else: kmer = 20

    if not pairing:
        input_file=intermediate_path + '{}/{}_nodonor.fastq'.format(input_name,input_name) # trimmed file name
        if not params['self-target_mask'][0] and metainfo['Genome']=='pace_ecoli' and metainfo['System']=='CAST': output_unmatched = params['Prefix']+'/alignment_input/' + input_name + '_input.fastq'
        else: output_unmatched = intermediate_path + '{}/{}_nospikein.fastq'.format(input_name,input_name)
        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched,
            'k={}'.format(kmer), 'hdist='+params['mismatch_tolerance'], 'literal=' + spikein_seq, 'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        retained_reads = read_len(output_unmatched)
        input_reads = read_len(input_file)

    elif pairing:
        input_file=intermediate_path + '{}/{}_nodonor'.format(input_name,input_name)
        if not params['self-target_mask'][0] and metainfo['Genome']=='pace_ecoli' and metainfo['System']=='CAST': output_unmatched = params['Prefix']+'/alignment_input/' + input_name + '_input' # Output file name for reads that don't have transposon end
        else: output_unmatched = intermediate_path + '{}/{}_nospikein'.format(input_name,input_name)
        
        bbduk_process = bbduk_func_paired(input_file,output_unmatched,spikein_seq)

        retained_reads = read_len('{}_1.fastq'.format(output_unmatched))
        input_reads = read_len("{}_1.fastq".format(input_file))
    
    spikein_reads = input_reads - retained_reads

    yields = get_yield(stats[1],spikein_reads)

    with open(intermediate_path + '{}/{}_nospikein_bbduk.txt'.format(input_name,input_name),'w') as f: f.write(bbduk_process.stderr)
    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [spikein_reads,0.000,retained_reads]
    return [spikein_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3],retained_reads]

def minicircle_finder(metainfo,stats,pairing=params['is_paired_end']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform Minicircle trimming for!\n")
        return ['NA', 'NA',0]
    input_name = metainfo['Sample']
    minicircle_seq = metainfo['Minicircle']
    if len(minicircle_seq)<20: kmer = len(minicircle_seq)
    else: kmer = 20


    if not pairing:
        input_file=intermediate_path + '{}/{}_nospikein.fastq'.format(input_name,input_name)
        output_unmatched = '{}/alignment_input/{}_input.fastq'.format(params['Prefix'],input_name)
        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched,
            'k={}'.format(kmer), 'hdist='+params['mismatch_tolerance'], 'literal=' + minicircle_seq, 'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        retained_reads = read_len(output_unmatched)
        input_reads = read_len(input_file)

    elif pairing:
        input_file=intermediate_path + '{}/{}_nospikein'.format(input_name,input_name)
        output_unmatched = '{}/alignment_input/{}_input'.format(params['Prefix'],input_name)

        bbduk_process = bbduk_func_paired(input_file,output_unmatched,minicircle_seq)

        retained_reads = read_len('{}_1.fastq'.format(output_unmatched))
        input_reads = read_len('{}_1.fastq'.format(input_file))
    
    minicircle_reads = input_reads - retained_reads

    with open(intermediate_path + '{}/{}_nominicircle_bbduk.txt'.format(input_name,input_name),'w') as f: f.write(bbduk_process.stderr)

    if retained_reads == 0: metainfo["Genome_mapping_reads"] = False
    return ['NA', 'NA',minicircle_reads]

def plasmid_finder(metainfo,stats,pairing=params['is_paired_end']):
    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform plasmid integration filtering for!\n")
        return ['NA', 'NA',0]
    input_name = metainfo['Sample']
    kmer = 20
    if not pairing:
        input_file=intermediate_path + '{}/{}_noself.fastq'.format(input_name,input_name)
        output_unmatched = '{}/alignment_input/{}_input.fastq'.format(params['Prefix'],input_name) # Output file name for reads that don't have transposon end
        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched, 'outm='+output_matched,
        'k={}'.format(kmer), 'hdist='+params['mismatch_tolerance'], 'ref='+index_folder+index_genomes['Transfected_plasmids'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        input_reads = read_len(input_file)
        retained_reads = read_len(output_unmatched)

    elif pairing:
        input_file=intermediate_path + '{}/{}_noself'.format(input_name,input_name)
        output_matched = '{}/alignment_input/{}_plasmid_reads'.format(params['Prefix'],input_name)
        output_unmatched = '{}/alignment_input/{}_input'.format(params['Prefix'],input_name) # Output file name for reads that don't have transposon end
        bbduk_process = subprocess.run(['bbduk.sh', 'in1={}_1.fastq in2={}_2.fastq'.format(input_file,input_file),
            'out1={}_1.fastq out2={}_2.fastq'.format(output_unmatched,output_unmatched),
            'outm1={}_1.fastq outm2={}_2.fastq'.format(output_matched, output_matched),
            'k={}'.format(kmer), 'hdist=2', 'ref='+index_folder+index_genomes['Transfected_plasmids'], '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)

        retained_reads = read_len('{}_1.fastq'.format(output_unmatched))
        input_reads = read_len("{}_1.fastq".format(input_file))
    
    donorintegration = input_reads - retained_reads

    with open(intermediate_path + '/{}/{}_plasmidfilters_bbduk.txt'.format(input_name,input_name),'w') as f: f.write(bbduk_process.stderr)

    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [donorintegration,0]
    return [donorintegration,retained_reads]

def self_targets(metainfo,stats,pairing=params['is_paired_end']): #NEED TO UPDATE 

    if metainfo["Genome_mapping_reads"] == False:
        print("No reads to perform Self-target trimming for!\n")
        return ['NA','NA',0]
    input_name = metainfo['Sample']
    crispr = metainfo['CRISPRarray']
    if not pairing:
        input_file=intermediate_path + '{}/{}_nospikein.fastq'.format(input_name,input_name)
        output_unmatched = intermediate_path + '{}/{}_noself.fastq'.format(input_name,input_name) # Output file name for reads that don't have transposon end
        bbduk_process = subprocess.run(['bbduk.sh', 'in=' + input_file, 'out=' + output_unmatched,
            'k=20', 'hdist='+params['mismatch_tolerance'],
            'literal=' + crispr, 'rcomp=f', '-Xmx{}g'.format(available_memory)], stderr=subprocess.PIPE, text=True)
        retained_reads = read_len(output_unmatched)
        input_reads = read_len(input_file)

    elif pairing:
        input_file=intermediate_path + '{}/{}_nospikein'.format(input_name,input_name) # trimmed file name
        output_unmatched = intermediate_path + '{}/{}_noself'.format(input_name,input_name) # Output file name for reads that don't have transposon end
        
        bbduk_process = bbduk_func_paired(input_file,output_unmatched,crispr)

        retained_reads = read_len('{}_1.fastq'.format(output_unmatched))
        input_reads = read_len("{}_1.fastq".format(input_file))

    with open(intermediate_path+'{}/{}_noself_bbduk.txt'.format(input_name,input_name),'w') as f: f.write(bbduk_process.stderr)

    removed_reads = input_reads - retained_reads
    yields = get_yield(stats[1],removed_reads)

    if retained_reads == 0:
        metainfo["Genome_mapping_reads"] = False
        return [input_reads,retained_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3]]
    return [removed_reads, yields.split(".")[0]+"."+yields.split(".")[1][:3], retained_reads]