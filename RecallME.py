#!/usr/bin/python #set the path for your env interpreter
#coding=utf-8

'''
Authors: Gianluca Vozza, Valentina Favalli, Emanuele Bonetti
v. 1.0
'''
# v. 1.0

# RecallME is free non-commercial software. 
# Users need to obtain the ANNOVAR licence by themselves. 
# Contact the Authors for commercial use.

#modules
import subprocess as sp
import sys
import os
import argparse
import tempfile
from os.path import isfile, join
#from pyfiglet import Figlet

#create ad hoc function for raising error if there's not vcf folder path
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)



#the parser  
parser = argparse.ArgumentParser(description='RecallMe is a software for computing variant calling pipelines metrics')
parser.add_argument('-q','--query_vcf', type=str,
                    help='query VCF file (or AVinput)', required=True)
parser.add_argument('--high_conf_bed', type=str,
                    help='Bed file of high confidence regions', required=False)
parser.add_argument('--query_format', type=str,
                    help='query file format (VCF or AVinput)', required=True)
parser.add_argument('--gt_format', type=str,
                    help='ground truth file format (VCF or AVinput)', required=True)
parser.add_argument('-f','--fasta', type=str,
                    help='fasta reference', required=True)
parser.add_argument('--genome_bed', type=str,
                    help='bedtools genome bed file with chromosome sizes ', required=True)
parser.add_argument('-b','--bam', type=str,
                    help='bam file of query', required=True)
parser.add_argument('-g','--ground_truth', type=str,
                    help='ground truth VCF file (or AVinput)', required=True)
parser.add_argument('-a','--annovar_dir', type=dir_path,
                    help='The path to the annovar directory', required=True)
parser.add_argument('-o','--out_dir', type=dir_path,
                    help='the path to output directory', required=True)
parser.add_argument('--vaf_query', type=float,
                    help='Set VAF threshold for query VCF (not required)', required=False)
parser.add_argument('--vaf_gt', type=float,
                    help='Set VAF threshold for ground truth VCF (not required)', required=False)
parser.add_argument('--caller', type=str,
                    help='Caller which produced the query VCF (GATK, TVC, Deepvariant, VarScan, LoFreq)', required=True)                   
args = parser.parse_args()

'''
def print_cool(text):
    cool_text = Figlet(font = "slant")
    os.system("cls")
    os.system('mode con: cols=75 lines=30')
    return str(cool_text.renderText(text))
    print(print_cool("RecallMe - © 2020"))
'''
print('Starting...')
print('Creating folder... Please wait')

#scripts folder
script_folder =  'scripts/'

#annovar folder
annovar_folder = args.annovar_dir

#Create AVinput directory
'''
command = 'mkdir ' + args.out_dir + 'AVinput'
process = sp.Popen(command, shell = True)
process.wait()
'''

#Create metrics directory
command = 'mkdir ' + args.out_dir + 'Metrics'
process = sp.Popen(command, shell = True)
process.wait()
print("Converting VCF files to AVinput files...")

if args.query_format == 'VCF':
    #Convert to AVinput query VCF
    command = 'perl ' + annovar_folder +  'convert2annovar.pl --format vcf4 ' + args.query_vcf + ' --includeinfo --outfile ' +  args.query_vcf.split('.')[0] + '.avinput'
    process = sp.Popen(command, shell = True)
    process.wait()
else:
    print('Query input is AVinput...Skipping conversion')

if args.gt_format == 'VCF':
        #Convert to AVinput GT VCF
        command = 'perl ' + annovar_folder + 'convert2annovar.pl --format vcf4 ' + args.ground_truth + ' --includeinfo --outfile ' + args.ground_truth.split('.')[0] + '.avinput' 
        process = sp.Popen(command, shell = True)
        process.wait()
else:
    print('GT input is AVinput...Skipping conversion')

print('Preparing for analysis...')


#Pass files to recalleR.r script
if args.vaf_query == None and args.vaf_gt == None:
    command = 'Rscript ' + script_folder + 'recaller.R -v ' + args.query_vcf.split('.')[0] + '.avinput' + ' -g '+ args.ground_truth.split('.')[0] + '.avinput' + ' --out ' + args.out_dir + 'Metrics/' + ' --caller ' + args.caller
    process = sp.Popen(command, shell = True)
    process.wait()
else:
    command = 'Rscript ' + script_folder + 'recaller.R -v ' + args.query_vcf.split('.')[0] + '.avinput' + ' -g ' + args.ground_truth.split('.')[0] + '.avinput'+ ' --query_vaf ' + args.vaf_query + ' --gt_vaf ' + args.vaf_gt + ' --out ' + args.out_dir + 'Metrics/' + ' --caller ' + args.caller
    process = sp.Popen(command, shell = True)
    process.wait()

print('Metrics generated.')
print('Starting mpileup...')

#FNs for mpileup
if os.path.isfile(args.out_dir +  'Metrics/' + 'FNs.txt'):
    command = 'awk \'{print $1,$2,$3}\' ' + args.out_dir + 'Metrics/' +  'FNs.txt > ' + args.out_dir + 'Metrics/' +  'FNs_for_pileup.txt'
    process = sp.Popen(command, shell = True)
    process.wait()

#mpileup
    command = 'bcftools mpileup -Ou -T ' + args.out_dir + 'Metrics/FNs_for_pileup.txt' + ' --gvcf 0 -f' + args.fasta + ' ' + args.bam + ' | bcftools call -m -Ov -o ' + args.out_dir + 'query_pileup.vcf'
    process = sp.Popen(command, shell=True)
    process.wait()

    print("Pileup ready!")
    print("Analyzing pileup...")
    print('Converting pileup into AVinput...')

    command = 'perl ' + annovar_folder + 'convert2annovar.pl --format vcf4 ' + args.out_dir + 'query_pileup.vcf --outfile ' + args.out_dir + 'query_pileup.avinput --includeinfo'
    process = sp.Popen(command, shell=True)
    process.wait()
    print('Starting analysis on actual False Negatives...')
    print('Loading into R...')
    print('Starting analysis...')

    if args.high_conf_bed == None:
        print('High Confidence Regions Bed file was not provided')
        print('Skipping Specificity computation...')

        command = 'Rscript ' + script_folder + 'pileup_recaller.R ' + '--pileup_file ' + args.out_dir + 'query_pileup.avinput ' + '--metrics ' + args.out_dir + 'metrics.txt ' +  '--TPs_table ' + args.out_dir + 'TPs.txt ' + '--FPs_table ' + args.out_dir + 'FPs.txt ' + '--FNs_table ' + args.out_dir + 'FNs.txt --out ' + args.out_dir
        process = sp.Popen(command, shell=True)
        process.wait()
        print('Tables generated!')

    else:
        print('High Confidence Regions Bed file was provided')
        print('Specificity computation...')
        command = 'bedtools genomecov -i ' + args.high_conf_bed + ' -g '+ args.genome_bed +  ' > ' + args.out_dir + 'genome_cov.txt'
        process = sp.Popen(command, shell=True)
        process.wait()

        #This is needed in order to create a bed file for subtraction by the high conf bed
        command = 'cat ' + args.out_dir + 'Metrics/TPs.txt ' + args.out_dir + 'Metrics/FNs.txt > ' + args.out_dir + 'bases_to_remove.txt'
        process = sp.Popen(command, shell=True)
        process.wait()
        
        command = 'grep -w genome ' + args.out_dir + 'genome_cov.txt | awk \'{print $3}\'  | awk \'FNR == 2 {print}\' > ' + args.out_dir + 'Metrics/bases.txt' 
        process = sp.Popen(command, shell=True)
        process.wait()

        command = 'Rscript ' + script_folder + 'pileup_recaller.R ' + '--pileup_file ' + args.out_dir + 'query_pileup.avinput ' + '--metrics ' + args.out_dir + 'Metrics/metrics.txt ' + '--bases ' + args.out_dir + 'Metrics/bases.txt ' + '--TPs_table ' + args.out_dir + 'Metrics/TPs.txt ' + '--FPs_table ' + args.out_dir + 'Metrics/FPs.txt ' + '--FNs_table ' + args.out_dir + 'Metrics/FNs.txt --out ' + args.out_dir + '/Metrics/'
        process = sp.Popen(command, shell=True)
        process.wait()

        print('Tables generated!')
        print('Removing temporary files')
        print('Please wait...')
        command = 'rm ' + args.out_dir + 'genome_cov.txt'
        process = sp.Popen(command, shell=True)
        process.wait()

        command = 'rm ' + args.out_dir + 'bases_to_remove.txt'
        process = sp.Popen(command, shell=True)
        process.wait()

        command = 'rm ' + args.out_dir + 'Metrics/bases.txt'
        process = sp.Popen(command, shell=True)
        process.wait()
    '''
        command = 'rm ' + args.out_dir + 'Metrics/TPs_tmp.txt'
        process = sp.Popen(command, shell=True)
        process.wait()
    '''

    command = 'rm ' + args.out_dir + 'Metrics/FNs_for_pileup.txt'
    process = sp.Popen(command, shell=True)
    process.wait()
    print("Analysis completed!")
    print("Thank you for using RecallMe - v.1.0")

else:
    print("Any FNs were generated")
    print("Metrics were not updated")
    print("Analysis completed!")
    print("Thank you for using RecallMe - v.1.0")
