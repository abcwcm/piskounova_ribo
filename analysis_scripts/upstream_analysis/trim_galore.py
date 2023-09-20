import os
import sys
import argparse
import logging
import re

'''Usage:trim_galore.py -f <path containing Raw files > -c <config file with software information> -a <single or pair end SE/PE>'''
#python trim_galore.py -f fastq.list -a 'SE/PE'
#trimm and output fastqc checking

def msg():
    return ''' This script runs adaptor filtering using trim-galore software.
               Exiting.... '''

def parse_arguments():
    '''Adding the command line arguments to be provided while running the script'''
    parser = argparse.ArgumentParser(description='Process command line arguments.', usage=msg())
    parser.add_argument('-f', '--file_list', metavar='FILE', required=True, help='fastq file list')
    parser.add_argument('-c', '--config', metavar='FILE', required=True, help='config file')
    parser.add_argument('-a', '--aim', metavar='SE/PE', required=True, help='SE/PE')
    args = parser.parse_args()
    return args

cmd_line_args = parse_arguments()
file_list = cmd_line_args.file_list
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




def trim_galore():

    my_config_dict = make_config_hash()

    sample_dictionary = make_sample_dict()
        

    for sample_name, dir_name in sample_dictionary.items():
        adaptor_option = my_config_dict['adaptor_option']
        sample_name2 = re.split("_S\d+", sample_name)[0]
        sample_name3 = re.split("_R\d+", sample_name)[0]
        if aim == 'SE':
            
            cmd = " echo trim_galore --fastqc " + adaptor_option + " -o " + dir_name + "/trimmed_out/ " + " --basename " + sample_name2 + " " + dir_name + "/" +sample_name
            return_val = os.system(cmd)

        if aim == 'PE':
            cmd = " echo trim_galore --paired --fastqc_args " + adaptor_option + " -o " + dir_name + "/trimmed_out/"  + " --retain_unpaired -r1 18 -r2 18" + " --basename " + sample_name2 + " " + dir_name + "/" +sample_name3 + "_R1_001.fastq.gz " +dir_name + "/" +sample_name3 + "_R2_001.fastq.gz " 
            return_val = os.system(cmd)

        
trim_galore()
