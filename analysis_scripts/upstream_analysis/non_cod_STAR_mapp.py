import os
import sys
import argparse
import logging
import re


'''Usage:STAR_mapping.py -f <path containing samples > -g <whole path containing genome> -c <config file with software information> -a <single or pair end SE/PE>'''
#python STAR_mapping.py -c config_mapping.yml -a PE -f sample.list -g /share/data/Illumina_Databases/human/star/genocode_GRCh38-v41/

#(needs to be the output of trimmomatic)
#ls -d /share/data/RNA_Seq/Matt_Cummings/Working/Trimmed_reads/output_* > sample.list
#cat sample.lst
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP22-PRE-Control-Null
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP30-PRE-Case-Null
#/share/data/RNA_Seq/ETT/star_Aug29_2022/test_raw_data/output_MP45-PRE-Case-Null

def msg():
    return ''' This script runs star mapping program to generate mapping alignemnts.
               Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-f', '--file_list', metavar='FILE', required=True, help='sample file list')
    parser.add_argument('-c', '--config', metavar='FILE', required=True, help='config file')
    parser.add_argument('-g', '--genome_path', metavar='DIR', required=True, help='genome folder whole path')
    parser.add_argument('-a', '--aim', metavar='SE/PE', required=True, help='SE/PE')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
file_list = cmd_line_args.file_list
genome_path = cmd_line_args.genome_path
aim = cmd_line_args.aim

def make_config_hash():
    config_dict = {}
    config_file = cmd_line_args.config
    config = open( config_file, 'r')

    for line in config:
        this_line = line.strip()
        this_key = this_line.split(":")[0]
        this_val = this_line.split(":")[1]
        config_dict[this_key] = this_val

    return config_dict


def make_sample_dict():

    sample_dict = {}

    #Iterate through the file list
    fastq_list = open(file_list, 'r')

    for fastq in fastq_list:
        this_fastq = fastq.strip()
        dir_name = os.path.dirname(this_fastq)
        file_name = os.path.basename(this_fastq)
        sample_name = file_name
        sample_dict[sample_name] = dir_name

    return sample_dict



def STAR_mapping():

    my_config_dict = make_config_hash()
    all_path_genome = os.path.abspath(genome_path)
    genome_dir = os.path.dirname(all_path_genome)
    
    genome_annot = all_path_genome.split(".")[0]
    #genome_annot = all_path_genome.split(".")[0]
    sample_dictionary = make_sample_dict()

    for sample_file, dir_name in sample_dictionary.items():
        num_cores = my_config_dict['num_cores']
        MultimappingMax = my_config_dict['max_multiMapping']
        EndsType = my_config_dict['EndsType']
        Unmapped_option = my_config_dict['Unmapped_option']
        alignIntronMax = my_config_dict['alignIntronMax']
        quantMode = my_config_dict['quantMode']
        sample_name2 = re.split("_val", sample_file)[0]
        #sample_name3= re.split("_val", sample_name2)[0]
        genome_dir = os.path.dirname(genome_path)


        if aim == 'SE':
            cmd = "echo  STAR --runThreadN " + num_cores + " --genomeDir " + genome_dir + \
            " --readFilesIn " + dir_name + "/" +sample_name2 +" --readFilesCommand zcat "  +\
            " --outFileNamePrefix " + dir_name + "/ouput/non_coding_" + sample_name2 +" --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax " + \
            MultimappingMax + " --quantMode " + quantMode +" --outReadsUnmapped " + Unmapped_option + " --alignIntronMax " + alignIntronMax +\
            " --limitBAMsortRAM 3726878714 --alignEndsType " + EndsType
            


        if aim == 'PE':

            cmd = "echo  STAR --runThreadN " + num_cores + " --genomeDir " + genome_dir + \
            " --readFilesIn " + dir_name + "/" +sample_name2 +"_val_1_outx.fq.gz " + dir_name + "/" +sample_name2 +"_val_2_outx.fq.gz " +" --readFilesCommand zcat "  +\
            " --outFileNamePrefix " + dir_name + "/ouput/non_coding_" + sample_name2 +" --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax " + \
            MultimappingMax + " --quantMode " + quantMode +" --outReadsUnmapped " + Unmapped_option + " --alignIntronMax " + alignIntronMax +\
            " --limitBAMsortRAM 1818334857 --alignEndsType " + EndsType
            
        return_val = os.system(cmd)

STAR_mapping()
