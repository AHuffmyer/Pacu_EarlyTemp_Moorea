# Analysis of Early Life History P. acuta Temp seqs   

### Ariana Huffmyer, 20220218  


# 1. Prepare data files on Andromeda  

Sym link files to working directory on Andromeda.  

```
cd /data/putnamlab/ashuffmyer
mkdir pacuta-tagseq
cd pacuta-tagseq
``` 

Sym link files into working directory.  
```
ln -s /data/putnamlab/shared/2022_Pacu_Moorea_TagSeq/*.fastq.gz /data/putnamlab/ashuffmyer/pacuta-tagseq/
```

This will create a symlink (blue) for each file in the working directory.  

Create a conda environment.  
```
conda create -n pacutatag
conda activate pacutatag
```

# 2. Generate MultiQC of raw files   

Run fastqc on the raw files and generate a multiQC report.  

```
mkdir fastqc_results 
cd ../ 
```

``` 
nano multiqc.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="multiqc_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="multiqc_output_script" #once your job is completed, any final job report comments will be put in this file

source /usr/share/Modules/init/sh # load the module function

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

for file in ./*fastq.gz
do
fastqc $file --outdir ./fastqc_results         
done

multiqc --interactive ./fastqc_results  

mv multiqc_report.html ./fastqc_results/pacuta_tagseq_raw_qc_multiqc_report.html #renames file
```

```
sbatch multiqc.sh
``` 

Output file to local computer and view (outside of Andromeda).  
```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/pacuta-tagseq/fastqc_results/pacuta_tagseq_raw_qc_multiqc_report.html ~/MyProjects/Pacu_EarlyTemp_Moorea/Pacu2021/Output/TagSeq/
```

# 3. Merge files for each sample across lanes.  


Merge across lanes.  

```
cd sequences 

nano merge.sh
``` 

```
#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab 

ls -1 Pact*R1*.gz | awk -F '_' '{print $1"_"$2}' | sort | uniq > ID

for i in `cat ./ID`; 

do cat "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz > "$i"_R1_cat.fastq.gz; 

done;

```

```
sbatch merge.sh
``` 

# 4. Conduct QC and filtering of sequences  


```
mkdir scripts
nano qc.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="qc_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="qc_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim
array1=($(ls *cat.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 
# fastqc the cleaned reads
        fastqc clean.${i}
done 

echo "Read trimming of adapters complete." $(date)

# Quality Assessment of Trimmed Reads

multiqc --interactive clean* #Compile MultiQC report from FastQC files - need to rename output file to "clean"

echo "Cleaned MultiQC report generated." $(date)

```

Option to do without the fastQC step: 
```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="qc_script_error" #if your job fails, the error report will be put in this file
#SBATCH --output="qc_output_script" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences

# load modules needed
module load fastp/0.19.7-foss-2018b
module load FastQC/0.11.8-Java-1.8
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# Make an array of sequences to trim
array1=($(ls *cat.fastq.gz)) 

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA) 
for i in ${array1[@]}; do
	fastp --in1 ${i} --out1 clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50 

done 

echo "Read trimming of adapters complete." $(date)
```

Rename and output clean multiQC file to local computer and view (outside of Andromeda).  

```
cd sequences 

mv multiqc_report.html ../fastqc_results/pacuta_tagseq_clean_qc_multiqc_report.html #renames file

scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/pacuta-tagseq/fastqc_results/pacuta_tagseq_clean_qc_multiqc_report.html ~/MyProjects/Pacu_EarlyTemp_Moorea/Pacu2021/Output/TagSeq/
```

# 5. HISAT2

Obtain reference genome assembly and gff annotation file.  

```
cd sequences 

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv1.genes.gff3.gz

wget http://cyanophora.rutgers.edu/Pocillopora_acuta/Pocillopora_acuta_HIv1.assembly.fasta.gz

```

Unzip gff file 

```
gunzip Pocillopora_acuta_HIv1.genes.gff3.gz
```

Create a script  

```
cd scripts

nano align.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="align_error" #if your job fails, the error report will be put in this file
#SBATCH --output="align_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences

# load modules needed
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.9-foss-2018b

#unzip reference genome
gunzip Pocillopora_acuta_HIv1.assembly.fasta.gz

# index the reference genome for Pocillopora acuta output index to working directory
hisat2-build -f /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences/Pocillopora_acuta_HIv1.assembly.fasta ./Pacuta_ref # called the reference genome (scaffolds)
echo "Reference genome indexed. Starting alingment" $(date)

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls clean*)) # call the clean sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	hisat2 -p 8 --dta -x Pacuta_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```

```
sbatch align.sh
```

This creates a .bam for every file. 

In the align_error file you can see the alignment stats - we had about 60-65% alignment in all samples.  

# 6. StringTie2

```
cd scripts
nano stringtie.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="stringtie_error" #if your job fails, the error report will be put in this file
#SBATCH --output="stringtie_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences

#load packages
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie

array=($(ls *.bam)) #Make an array of sequences to assemble
 
for i in ${array[@]}; do #Running with the -e option to compare output to exclude novel genes. Also output a file with the gene abundances
        sample_name=`echo $i| awk -F [_] '{print $1"_"$2"_"$3}'`
	stringtie -p 8 -e -B -G Pocillopora_acuta_HIv1.genes.gff3 -A ${sample_name}.gene_abund.tab -o ${sample_name}.gtf ${i}
        echo "StringTie assembly for seq file ${i}" $(date)
done
echo "StringTie assembly COMPLETE, starting assembly analysis" $(date)

```

-p means number of threads/CPUs to use (8 here)

-e means only estimate abundance of given reference transcripts (only genes from the genome) - dont use if using splice variance aware to merge novel and ref based.  

-B means enable output of ballgown table files to be created in same output as GTF

-G means genome reference to be included in the merging 

```
sbatch stringtie.sh
```

This will make a .gtf file for each sample. 

# 7. Prep DE

Grab the python script from another AH project. 

```
cp /data/putnamlab/ashuffmyer/pairs-rnaseq/prepDE.py /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences
```

```
cd scripts
nano prepDE.sh
```

```
#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=ashuffmyer@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="prepDE_error" #if your job fails, the error report will be put in this file
#SBATCH --output="prepDE_output" #once your job is completed, any final job report comments will be put in this file
#SBATCH -D /data/putnamlab/ashuffmyer/pacuta-tagseq/sequences

#load packages
module load Python/2.7.15-foss-2018b #Python
module load StringTie/2.1.4-GCC-9.3.0

#Transcript assembly: StringTie
module load GffCompare/0.12.1-GCCcore-8.3.0

#Transcript assembly QC: GFFCompare

#make gtf_list.txt file
ls Pact*.gtf > gtf_list.txt

stringtie --merge -p 8 -G Pocillopora_acuta_HIv1.genes.gff3 -o Pacuta_merged.gtf gtf_list.txt #Merge GTFs to form $
echo "Stringtie merge complete" $(date)

gffcompare -r Pocillopora_acuta_HIv1.genes.gff3 -G -o merged Pacuta_merged.gtf #Compute the accuracy and pre$
echo "GFFcompare complete, Starting gene count matrix assembly..." $(date)

#make gtf list text file
for filename in Pact*.gtf; do echo $filename $PWD/$filename; done > listGTF.txt

python prepDE.py -g Pacuta_gene_count_matrix.csv -i listGTF.txt #Compile the gene count matrix
echo "Gene count matrix compiled." $(date)
```

```
sbatch prepDE.sh
```


Move gene count matrix off of server  

```
scp ashuffmyer@ssh3.hac.uri.edu:/data/putnamlab/ashuffmyer/pacuta-tagseq/sequences/Pacuta_gene_count_matrix.csv ~/MyProjects/Pacu_EarlyTemp_Moorea/Pacu2021/Output/TagSeq/
```


