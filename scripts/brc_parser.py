import os
import argparse
from collections import defaultdict
import sys
import csv
import re
import pandas as pd



def main():

		'''
				 Given bam-readcount output
					produces output is a set of comma
					delimited columns, indicating indels and snvs
					strand counts:w etc..
		'''
		parser = argparse.ArgumentParser(
				formatter_class=argparse.RawDescriptionHelpFormatter,
				description='bam-read count parser')
		parser.add_argument("infile",
						type=str,
							action='store',
							help='inputfile from bam count ouput')

		args = parser.parse_args()

		parser_readcount(args)




def parser_readcount(args):
	ids = defaultdict(list) # initiate a default dict(list)

	with open(args.infile, 'r') as f:
			reader = csv.reader(f, delimiter='\t')
			for row in reader:
					chrm, key, ref, dp = row[:4]
					#print chrm, key, ref, dp
					for elems in row[4:]:
							elems = elems.split(':')
							#print elems

							for value in  elems[1]:
									if elems[1] != '0':
										pb = round(float(elems[5]) / float(elems[1]), 2)
										vaf = round(float(elems[1]) / int(dp),2)
									else:
										pb = '0'
										vaf = '0'
							_bas = (elems[0])
							#print(_bas)
							#for _base in _bas:
									#print _base
							if _bas.startswith("+") == True and _bas != ref :
												#print ('inertion:' , _base)
									_mut = 'ins'
							elif _bas.startswith("-") == True and _bas != ref:
									_mut = 'del'
							elif _bas != ref:
										_mut = 'snv'
							elif _bas in ref:
									_mut = 'no-mutation'
							else:
									_mut = 'NA'

							if elems[0] != '=':
											ids[key].append({
										'chr': chrm,
										'ref':ref,
										'depth' : dp,
										'base': '(' + elems[0]+ ')',
										'count': elems[1],
										'positive_strand': elems[5],
										'negative_strand': elems[6],
										'percent_bias': pb,
										'vaf': vaf,
										'mutation' : _mut

										})


	with open(args.infile[:-4]+'_parsed.csv', 'w+') as f:
				writer = csv.writer(f, delimiter=",")
				keys = ["chr","ref",
								"depth", "count",
								"base", "positive_strand", "negative_strand",
								"percent_bias",
								"vaf", 'mutation']
				writer.writerow(["position"] + keys)
				for k, vl in ids.items():
						for v in vl:
							writer.writerow([k] + [v[key] for key in keys])


if __name__ == '__main__':
		main()
