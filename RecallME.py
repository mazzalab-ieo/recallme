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

#Function that raises an error if the string passed to an argument is not a dir
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)



#Parsing arguments  
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
                    help='Set VAF threshold for query VCF (optional)', required=False)
parser.add_argument('--vaf_gt', type=float,
                    help='Set VAF threshold for ground truth VCF (optional)', required=False)
parser.add_argument('--caller', type=str,
                    help='Caller which produced the query VCF (GATK, TVC, Deepvariant, VarScan2, LoFreq, VarDict)', required=True)
parser.add_argument('--report', type=str,
                    help='Produce HTML report (yes or no, not required)', required=False)                  
args = parser.parse_args()

#Set folders
#scripts folder
script_folder =  'scripts/'

#annovar folder
annovar_folder = args.annovar_dir

#Create a folder for Metrics
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

print('Defining variant types...')

## variant types definition
command = 'python ' + script_folder + 'variant_type.py --query ' + args.query_vcf.split('.')[0] + '.avinput' + ' --gt ' + args.ground_truth.split('.')[0] + '.avinput'
process = sp.Popen(command, shell = True)
process.wait()

#Pass files to recalleR.r script
if args.vaf_query == None and args.vaf_gt == None:
    command = 'Rscript ' + script_folder + 'recaller.R -v ' + args.query_vcf.split('.')[0] + '_var_type.avinput' + ' -g '+ args.ground_truth.split('.')[0] + '_var_type.avinput' + ' --out ' + args.out_dir + 'Metrics/' + ' --caller ' + args.caller
    process = sp.Popen(command, shell = True)
    process.wait()
else:
    command = 'Rscript ' + script_folder + 'recaller.R -v ' + args.query_vcf.split('.')[0] + '_var_type.avinput' + ' -g ' + args.ground_truth.split('.')[0] + '_var_type.avinput' + ' --query_vaf ' + args.vaf_query + ' --gt_vaf ' + args.vaf_gt + ' --out ' + args.out_dir + 'Metrics/' + ' --caller ' + args.caller
    process = sp.Popen(command, shell = True)
    process.wait()

print('Metrics generated.')
print('Starting mpileup...')

#FNs for mpileup on snv
if os.path.isfile(args.out_dir +  'Metrics/' + 'FNs_snv.txt'):
    command = 'awk \'{print $1,$2,$3}\' ' + args.out_dir + 'Metrics/' +  'FNs_snv.txt > ' + args.out_dir + 'Metrics/' +  'FNs_snv_for_pileup.txt'
    process = sp.Popen(command, shell = True)
    process.wait()
    #mpileup for snv
    command = 'bcftools mpileup -Ou -T ' + args.out_dir + 'Metrics/FNs_snv_for_pileup.txt' + ' --gvcf 0 -f' + args.fasta + ' ' + args.bam + ' | bcftools call -m -Ov -o ' + args.out_dir + 'query_pileup_snv.vcf'
    process = sp.Popen(command, shell=True)
    process.wait()
    #annovar conversion SNVs
    command = 'perl ' + annovar_folder + 'convert2annovar.pl --format vcf4 ' + args.out_dir + 'query_pileup_snv.vcf --outfile ' + args.out_dir + 'query_pileup_snv.avinput --includeinfo'
    process = sp.Popen(command, shell=True)
    process.wait()
#FNs for mpileup on indels
if os.path.isfile(args.out_dir +  'Metrics/' + 'FNs_indel.txt'):
    command = 'awk \'{print $1,$2,$3}\' ' + args.out_dir + 'Metrics/' +  'FNs_indel.txt > ' + args.out_dir + 'Metrics/' +  'FNs_indel_for_pileup.txt'
    process = sp.Popen(command, shell = True)
    process.wait()
    #mpileup for indels
    command = 'bcftools mpileup -Ou -T ' + args.out_dir + 'Metrics/FNs_indel_for_pileup.txt' + ' --gvcf 0 -f' + args.fasta + ' ' + args.bam + ' | bcftools call -m -Ov -o ' + args.out_dir + 'query_pileup_indel.vcf'
    process = sp.Popen(command, shell=True)
    process.wait()

    print("Pileup ready!")
    print("Analyzing pileup...")
    print('Converting pileup into AVinput...')
    #annovar conversion INDELs
    command = 'perl ' + annovar_folder + 'convert2annovar.pl --format vcf4 ' + args.out_dir + 'query_pileup_indel.vcf --outfile ' + args.out_dir + 'query_pileup_indel.avinput --includeinfo'
    process = sp.Popen(command, shell=True)
    process.wait()

    print('Starting analysis on actual False Negatives...')
    print('Loading into R...')
    print('Starting analysis...')

if args.high_conf_bed == None:
    print('High Confidence Regions Bed file was not provided')
    print('Skipping Specificity computation...')
    command = 'Rscript ' + script_folder + 'pileup_recaller.R ' + '--pileup_file_snv ' + args.out_dir + 'query_pileup_snv.avinput ' + '--pileup_file_indel ' + args.out_dir + 'query_pileup_indel.avinput ' + '--metrics_snv ' + args.out_dir + 'Metrics/metrics_snv.txt ' + '--metrics_indel ' + args.out_dir + 'Metrics/metrics_indel.txt ' +  '--TPs_table_snv ' + args.out_dir + 'Metrics/TPs_snv.txt ' + '--TPs_table_indel ' + args.out_dir + 'Metrics/TPs_indel.txt ' + '--FPs_table_snv ' + args.out_dir + 'Metrics/FPs_snv.txt ' + '--FPs_table_indel ' + args.out_dir + 'Metrics/FPs_indel.txt ' + '--FNs_table_snv ' + args.out_dir + 'Metrics/FNs_snv.txt' + ' --FNs_table_indel ' + args.out_dir + 'Metrics/FNs_indel.txt' + ' --out ' + args.out_dir + 'Metrics/'
    process = sp.Popen(command, shell=True)
    process.wait()
    print('Tables generated!')
    command = 'Rscript ' + script_folder + 'pileup_recaller.R ' + '--pileup_file_snv ' + args.out_dir + 'query_pileup_snv.avinput ' + '--pileup_file_indel ' + args.out_dir + 'query_pileup_indel.avinput ' + '--metrics_snv ' + args.out_dir + 'Metrics/metrics_snv.txt ' + '--metrics_indel ' + args.out_dir + 'Metrics/metrics_indel.txt ' +  '--TPs_table_snv ' + args.out_dir + 'Metrics/TPs_snv.txt ' + '--TPs_table_indel ' + args.out_dir + 'Metrics/TPs_indel.txt ' + '--FPs_table_snv ' + args.out_dir + 'Metrics/FPs_snv.txt ' + '--FPs_table_indel ' + args.out_dir + 'Metrics/FPs_indel.txt ' + '--FNs_table_snv ' + args.out_dir + 'Metrics/FNs_snv.txt' + ' --FNs_table_indel ' + args.out_dir + 'Metrics/FNs_indel.txt' + ' --out ' + args.out_dir + 'Metrics/'
    process = sp.Popen(command, shell=True)
    process.wait()
    print('Tables generated!')
else:
    print('High Confidence Regions Bed file was provided')
    print('Specificity computation...')
    command = 'bedtools genomecov -i ' + args.high_conf_bed + ' -g '+ args.genome_bed +  ' > ' + args.out_dir + 'genome_cov.txt'
    process = sp.Popen(command, shell=True)
    process.wait()

    command = 'grep -w genome ' + args.out_dir + 'genome_cov.txt | awk \'{print $3}\'  | awk \'FNR == 2 {print}\' > ' + args.out_dir + 'Metrics/bases.txt' 
    process = sp.Popen(command, shell=True)
    process.wait()

    command = 'Rscript ' + script_folder + 'pileup_recaller.R ' + '--pileup_file_snv ' + args.out_dir + 'query_pileup_snv.avinput ' + '--pileup_file_indel ' + args.out_dir + 'query_pileup_indel.avinput ' + '--metrics_snv ' + args.out_dir + 'Metrics/metrics_snv.txt ' + '--metrics_indel ' + args.out_dir + 'Metrics/metrics_indel.txt ' + '--bases ' + args.out_dir + 'Metrics/bases.txt ' + '--TPs_table_snv ' + args.out_dir + 'Metrics/TPs_snv.txt ' + '--TPs_table_indel ' + args.out_dir + 'Metrics/TPs_indel.txt ' + '--FPs_table_snv ' + args.out_dir + 'Metrics/FPs_snv.txt ' + '--FPs_table_indel ' + args.out_dir + 'Metrics/FPs_indel.txt ' + '--FNs_table_snv ' + args.out_dir + 'Metrics/FNs_snv.txt' + ' --FNs_table_indel ' + args.out_dir + 'Metrics/FNs_indel.txt' + ' --out ' + args.out_dir + 'Metrics/'
    process = sp.Popen(command, shell=True)
    process.wait()
    '''
    print('Tables generated!')
    print('generating ReportME...')
    command = 'Rscript -e \"rmarkdown::render(' + script_folder + '\'ReportME.Rmd\', params = list(\'' + args.query_vcf.split('.')[0] + '_var_type.avinput' + '\',' + '\'' + args.out_dir + 'Metrics/TPs_snv.txt\',' + ' \'' + args.out_dir + 'Metrics/TPs_indel.txt\',' + ' \'' + args.out_dir + 'Metrics/FPs_snv.txt\','+ ' \'' + args.out_dir + 'Metrics/FPs_indel.txt\',' + ' \'' + args.out_dir + 'Metrics/FNs_snv.txt\',' + ' \'' + args.out_dir + 'Metrics/FNs_indel.txt\',' + ' \'' + args.out_dir + 'Metrics/metrics_snv_updated.txt\',' + ' \'' + args.out_dir + 'Metrics/metrics_indel_updated.txt\',' + '\''+ args.caller + '\')'
    process = sp.Popen(command, shell=True)
    process.wait()
    '''
    print('Tables generated!')
    print('generating ReportME...')
    if args.report == 'yes':
        command = 'Rscript ' + script_folder + 'report_gen.R  ' + args.query_vcf.split('.')[0] + '_var_type.avinput ' + args.out_dir + 'Metrics/TPs_snv.txt ' + args.out_dir + 'Metrics/TPs_indel.txt ' + args.out_dir + 'Metrics/FPs_snv.txt ' + args.out_dir + 'Metrics/FPs_indel.txt ' + args.out_dir + 'Metrics/FNs_snv.txt ' + args.out_dir + 'Metrics/FNs_indel.txt ' + args.out_dir + 'Metrics/metrics_snv_updated.txt ' + args.out_dir + 'Metrics/metrics_indel_updated.txt ' + args.caller
        process = sp.Popen(command, shell=True)
        process.wait()

        command = 'mv ' + script_folder + 'ReportME.html ' + args.out_dir + 'Metrics/'
        process = sp.Popen(command, shell=True)
        process.wait()
    else:
        print('Skipping report...')
print('Removing temporary files')
print('Please wait...')
command = 'rm ' + args.out_dir + 'genome_cov.txt'
process = sp.Popen(command, shell=True)
process.wait()

command = 'rm ' + args.out_dir + 'Metrics/bases.txt'
process = sp.Popen(command, shell=True)
process.wait()

command = 'rm ' + args.out_dir + 'Metrics/FNs_snv_for_pileup.txt'
process = sp.Popen(command, shell=True)
process.wait()
print("Analysis completed!")
print("Thank you for using RecallMe - v.1.0")
