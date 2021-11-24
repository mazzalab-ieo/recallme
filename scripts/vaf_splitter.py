#coding=utf-8
'''
Authors: G. Vozza, E. Bonetti
v. 0.1

This script retrieve VAF information from thermo fisher VCFs
'''
# modules
import os
import sys

out_file = open(sys.argv[2], 'w') #open new VCF file to write
with open(sys.argv[1], "r") as query: #open query VCF to split
    for line in query:
        if line.startswith('#'): #add VCF header to the new file
            out_file.write(line)
            continue
        split_line = line.strip().split('\t')
        info = split_line[7] 
        info_gt = split_line[9] # le informazioni su genotipo, AF ecc. sono contenute nella colonna numero 9
        alt = split_line[4]
        if ',' in alt:
            alt_array = alt.split(',') #create an array with all the alterations in the same position
            n_alternative = len(alt_array) #get the number of alterations in the position as an integer
            info_split = info.split(';')
            info_gt_split = info_gt.split(':')
            for i in range(0, n_alternative):
                string_to_write = '\t'.join(split_line[0:4]) + '\t' + alt_array[i] + '\t' + split_line[5] + '\t' + split_line[6] + '\t' #write many rows as the number of the alterations at the same position
                for element in info_split:
                    split_element = element.split('=')
                    info_symbol = split_element[0]
                    if ',' in split_element[1]:
                        if len(split_element[1].split(',')) < n_alternative: #if elements number is not the same of the number of alternatives then assign the whole string to all the alternatives
                            info_value = split_element[1]
                        else:
                            info_value = split_element[1].split(',')[i]
                        string_to_write = string_to_write + info_symbol + '=' + info_value + ';'
                    else:
                        info_value = split_element[1]
                        string_to_write = string_to_write + info_symbol + '=' + info_value + ';'
                string_to_write = string_to_write[:-1] #remove ';' at the end of the line
                string_to_write = string_to_write + '\t' + split_line[8] + '\t' + info_gt + '\n'
                out_file.write(string_to_write)
        else:
            out_file.write(line)
out_file.close()                

