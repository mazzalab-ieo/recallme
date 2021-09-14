#coding=utf-8
'''
Author: G. Vozza
v. 0.1
'''
#modules
import os
import argparse

parser = argparse.ArgumentParser(description='This program assign SNV or INDEL class to variants in AVinput for RecallME') #program description
parser.add_argument('--query', type=str,
                    help='query AVinput', required=True)# create the argument
parser.add_argument('--gt', type=str,
                    help='gt AVinput', required=True)# create the argument
args = parser.parse_args()

with open(args.query, 'r') as query:
    with open(args.query.split(".")[0] + '_var_type.avinput', 'w') as write_query:
        for line in query:
            total_line = line.strip('\n')
            ref = line.strip().split('\t')[3]
            alt = line.strip().split('\t')[4]
            if len(ref) == len(alt) and ref != "-" and alt != "-":
                write_query.write(total_line + '\t' + "SNV" + '\n') 
            elif len(ref) != len(alt) and ref == "-" or alt == "-":
                write_query.write(total_line + '\t' + "INDEL" + '\n')
            elif len(ref) == len(alt) and ref == "-" or alt == "-":
                write_query.write(total_line + '\t' + "INDEL" + '\n')
            else:
                write_query.write(total_line + '\t' + "OTHER" + '\n')
        write_query.close()

with open(args.gt, 'r') as gt:
    with open(args.gt.split(".")[0] + '_var_type.avinput', 'w') as write_gt:
        for line in gt:
            total_line = line.strip('\n')
            ref = line.strip().split('\t')[3]
            alt = line.strip().split('\t')[4]
            if len(ref) == len(alt) and ref != "-" and alt != "-":
                write_gt.write(total_line + '\t' + "SNV" + '\n') 
            elif len(ref) != len(alt) and ref == "-" or alt == "-":
                write_gt.write(total_line + '\t' + "INDEL" + '\n')
            elif len(ref) == len(alt) and ref == "-" or alt == "-":
                write_gt.write(total_line + '\t' + "INDEL" + '\n')
            else:
                write_gt.write(total_line + '\t' + "OTHER" + '\n')                
        write_gt.close()
        
