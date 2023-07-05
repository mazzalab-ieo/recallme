RecallME: a tool for VCF files comparison and efficient validation method
=============
![latest release](https://img.shields.io/github/v/release/mazzalab-ieo/recallme)
![Watch](https://img.shields.io/github/watchers/mazzalab-ieo/recallme?label=Watch)
![Stars](https://img.shields.io/github/stars/mazzalab-ieo/recallme?style=social)
![](https://img.shields.io/static/v1?label=Platform&message=Linux&color=lightgrey)
![](https://img.shields.io/static/v1?label=Dependencies&message=singularity,annovar&color=lightgrey)
![](https://img.shields.io/static/v1?label=Container&message=Docker&color=lightgrey)


RecallME is a containerized python tool for VCF file comparison and variant calling pipelines benchmarking and optimization.

### Quick start
To run RecallME you need to activate the conda environment and to install singularity (> v.3) and ANNOVAR (it requires a license free of charge for non-commercial purposes). For more information about ANNOVAR, please visit the <a href="https://annovar.openbioinformatics.org/en/latest/" target="_blank"> official website</a>


### Installation and setup
Once the repository has been cloned, the user has to create the conda environment by typing the following command:
```
conda env create --file env.txt --name {name_env}
```

### Basic commands
To run a comparison, type:
```
python RecallME.py -q $QUERY_VCF \
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
```
python RecallME.py --help
```
### Tutorial 
To run the example, open a terminal and type:
```
cd recallme/
bash run_example.sh
```

<code>run_example.sh</code> produces an outdir/ folder in the current directory as follows:<br>

<img src="www/diag.png" alt="diag" width="500" height="333">
              
After the output has been generated, go to the <a href="https://translational-oncology-lab.shinyapps.io/recallme/" target="_blank">web app</a> and load the AVinputs (query_split_var_type.avinput + ground_truth_var_type.avinput) and the rds object to visualize your results and to optimize your pipeline.

For additional information, please check the <a href="https://translational-oncology-lab.shinyapps.io/recallme/" target="_blank">about</a> page for interrogating the manual.

### License
**MIT License**

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
