#coding=utf-8
'''
Author: G. Vozza
v. 0.1
'''
#modules
import os
import argparse

parser = argparse.ArgumentParser(description='This program assign SNV or INDEL class to variants in AVinput for RecallME') #program description
parser.add_argument('--indel_avinput', type=str,
                    help='query AVinput', required=True)# create the argument
parser.add_argument('--converted_avinput', type=str,
                    help='converted AVinput', required=True)# create the argument
args = parser.parse_args()
with open(args.indel_avinput, 'r') as pileup:
    with open(args.converted_avinput, 'a') as new_file:
        for line in pileup:
            split_line = line.strip().split('\t')
            chr = split_line[0]
            start = split_line[1]
            end = split_line[2]
            ref = split_line[3]
            alt = split_line[4]
            chr_vcf = split_line[5]
            pos_vcf = split_line[6]
            id_vcf = split_line[7]
            ref_vcf = split_line[8]
            alt_vcf = split_line[9]
            info1 = split_line[10]
            info2 = split_line[11]
            info3 = split_line[12]
            info4 = split_line[13]
            info5 = split_line[14]
            if ref == '0':
                ref = '-'
            elif alt == '0':
                alt = '-'
            new_line = chr + '\t' + start + '\t' + end + '\t' + ref + '\t' + alt + '\t' + chr_vcf + '\t' + pos_vcf + '\t' + id_vcf + '\t' + ref_vcf + '\t' + alt_vcf + '\t' + info1 + '\t' + info2 + '\t' + info3 + '\t' + info4 + '\t' + info5 + '\n'
            new_file.write(new_line)

new_file.close()