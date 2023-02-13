## RecallME: a tool for VCF files comparison and efficient validation method
RecallME is a condarized python tool for VCF files comparison and variant calling pipelines benchmarking and optimization.

## Installation and setup
Once downloaded the repository, user has to install the conda environment from the yml file by typing the following command
```
conda env create -f RecallME_0.1.yml 
```
To run RecallME you need to activate the conda environment and to install ANNOVAR (it requires a license which is free for non-commercial use).
For more information about ANNOVAR, please visit https://annovar.openbioinformatics.org/en/latest/

## Basic commands
To run a comparison, type:
```
python -q $QUERY_VCF \
--high_conf_bed $BED \
--query_format VCF \
-g $GROUND_TRUTH \
--gt_format VCF \
-f $FASTA \
--genome_bed $GENOME_BED \
-b $BAM \
-a $ANNOVAR_DIR \
-o $OUT_DIR \
--caller $CALLER_NAME 
```

* **GENOME_BED** is a bedtools specific bed file with chromosome sizes in bps (i.e. chr1 n) or a .fai index.
* **BAM** file for query
* **ANNOVAR_DIR** the directory where is stored convert2annovar.pl which is required for running RecallME
* **CALLER_NAME** the name of the caller used to produce the query VCF file (GATK, TVC, Deepvariant, VarScan, LoFreq)

For more info on required and not required commands, please type:
python RecallME.py --help

## Pipeline optimization
Within conda env start RecallMEapp.Rmd shiny app and start to refine pipeline parameters thresholds.

## Info and comments
If you need help or you have comments and tips for improving RecallME, please send a mail to gianluca.vozza@ieo.it

## License

 RecallME is a free non-commercial software. Users need to obtain the ANNOVAR license by themselves.
 
 ## Reference
 
RecallME is currently under submission to a peer-reviewed journal.
