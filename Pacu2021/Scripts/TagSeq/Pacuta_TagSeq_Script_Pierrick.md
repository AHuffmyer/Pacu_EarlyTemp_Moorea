# Analysis of Early Life History P. acuta Temp seqs   

### Ariana Huffmyer, Pierrick Harnay 20220308 


# 1. Prepare data files on Andromeda  

Log onto Andromeda and go to the data folder 
 

```
cd /data/putnamlab/pharnay
cd Pacuta_Thermal/raw
ls 
```

Give Ariana editing permissions for your files.  

```
chmod u=rwx,g=rwx,o=rwx,a=rwx -R Pacuta_Thermal
```

In this folder, we have our fastq.gz files for all of our samples.  

We first need to look at the quality information for our sequences.  

# 2. Generate MultiQC of raw files   

Run fastqc on the raw files and generate a multiQC report.  

```
mkdir fastqc_results 
cd ../ 
```

Now we need to write a script.  

``` 
nano multiqc.sh
```

Copy in this information.  

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --mail-user=pierrick_harnay@uri.edu #your email to send notifications
#SBATCH --account=putnamlab                  
#SBATCH --error="multiqc_error_script" #if your job fails, the error report will be put in this file
#SBATCH --output="multiqc_output_script" #once your job is completed, any final job report comments will be put in this file

# load the necessary programs
source /usr/share/Modules/init/sh # load the module function

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-intel-2020a-Python-3.8.2

# run the fastQC on every sample
for file in ./raw/*fastq.gz
do
fastqc $file --outdir ./fastqc_results         
done

multiqc --interactive ./fastqc_results  

mv multiqc_report.html ./fastqc_results/pacuta_tagseq_raw_qc_multiqc_report.html #renames file
```

```
sbatch multiqc.sh
``` 

Check the status of the job.  

```
squeue -u pierrick_harnay
```

Once the job is done, we want to view the newly made file. 

Open a new terminal that is not logged into Andromeda. Output file to local computer and view (outside of Andromeda).  

```
scp pierrick_harnay@ssh3.hac.uri.edu:/data/putnamlab/pharnay/Pacuta_Thermal/fastqc_results/pacuta_tagseq_raw_qc_multiqc_report.html ~Desktop
```
