import sys
with open(sys.argv[1], 'r') as vcf:
    with open(sys.argv[2] + 'query_pileup_indel_dash.vcf', 'a') as new_vcf:
        for line in vcf:
            if line.startswith('#'):
                strip_line = line.strip()
                new_vcf.write(strip_line + '\n')
                continue
            else:
                split_line = line.strip().split('\t')
                chr = split_line[0]
                pos = split_line[1]
                id = split_line[2]
                ref = split_line[3]
                alt = split_line[4]
                qual = split_line[5]
                filter = split_line[6]
                info = split_line[7]
                format = split_line[8]
                sample = split_line[9]
                if '.' in ref:
                    ref = '-'
                elif '.' in alt:
                    alt = '-'
                new_line = chr + '\t' + pos + '\t' + id + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filter + '\t' + info + '\t' + format + '\t' + sample 
                new_vcf.write(new_line + '\n')
    new_vcf.close()