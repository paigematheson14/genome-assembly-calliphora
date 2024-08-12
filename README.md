# Genome assembly - Calliphora sp. 
This repository will detail the workflow I used to create reference genomes for four Calliphora species (C. quadrimaculata, C. hilli, C. stygia, and C. vicina) using PromethION nanopore data and Illumina data. 

# Folder organisation
Creating a genome is a complex process that generates many different files. Thus, it is important to keep folders organised in a way that makes sense. Meeran has a great github page that is dedicated to this and I recommend structuring folders in this manner - https://github.com/meeranhussain/QuickTips

# Sequenced data - PromethION data
My promethion data arrived to me in BLOW5 format. There were two libraries (pool A and pool B) which contained data for four samples (each sample representing one of the four aforementioned species). There were 96 barcodes (in blow5 format) in each library, however, we only really care about the first four barcodes (as these represent our samples and contain the bulk of the data we need). I think that the remaining 92 barcodes just contain junk but I'm not 100% sure on that fact. 

- NB01	GR-AC-CQ-50 (i.e., calliphora quadrimaculata) 
- NB02	GR-AC-CH-259 (i.e., calliphora hilli)
- NB03	GR-AC-CS-235 (i.e., calliphora stygia)
- NB04	GR-AC-CV-196 (i.e., calliphora vicina)

My Nanopore data came to me in a basecalled format (i.e., FASTQ). However, Nanopore currently uses Guppy to do basecalling which is not as effective as Dorado (Nanopore is going to switch to using Dorado eventually as it is actually better). Therefore, we decided that we would start with the raw BLOW5 files and do the basecalling step ourselves. 

# 1. Convert BLOW5/SLOW5 files to POD5

I used blue-crab (https://github.com/Psy-Fer/blue-crab) to make this conversion: 

```blue-crab s2p barcodes -d POD5_BARCODES_17JUNE```

# 2. Basecalling 

I used DORADO (https://github.com/nanoporetech/dorado). You first need to create a sample file (in .csv format) that contains metadata for your samples like so: 

| #  | kit                | experiment_id | flow_cell_id | barcode   | alias        |
|----|--------------------|---------------|--------------|-----------|--------------|
| 1  | SQK-MLK111-96-XL   | not_set       | PAK17764     | barcode01 | GR_AC_CQ_50  |
| 2  | SQK-MLK111-96-XL   | not_set       | PAK17764     | barcode02 | GR_AC_CH_259 |
| 3  | SQK-MLK111-96-XL   | not_set       | PAK17764     | barcode03 | GR_AC_CS_235 |
| 4  | SQK-MLK111-96-XL   | not_set       | PAK17764     | barcode04 | GR_AC_CV_196 |

For Pool A, and 

| #  | kit                | experiment_id | flow_cell_id | barcode   | alias        |
|----|--------------------|---------------|--------------|-----------|--------------|
| 1  | SQK-MLK111-96-XL   | not_set       | PAK33653     | barcode01 | GR_AC_CQ_50  |
| 2  | SQK-MLK111-96-XL   | not_set       | PAK33653     | barcode02 | GR_AC_CH_259 |
| 3  | SQK-MLK111-96-XL   | not_set       | PAK33653     | barcode03 | GR_AC_CS_235 |
| 4  | SQK-MLK111-96-XL   | not_set       | PAK33653     | barcode04 | GR_AC_CV_196 |

For Pool B. 

NOTE: you will need the kit name and experiment ID to match what is in your data. If you didn't recieve this information from your sequencer (like I), it is possible to retrieve this information from the POD5 files using the POD5 tool (https://github.com/nanoporetech/pod5-file-format). 

The CSV file needs to be in the same folder as the POD5 files. If you have two libraries to basecall, then you need to have two folders and run each library seperately. Here is the slurm script for basecalling (BTW you only need the '#SBATCH --gpus-per-node=A100:1' line for the basecalling part, so can remove from subsequent slurms):

```
#!/bin/bash -e

#SBATCH --account=uow03920
#SBATCH --job-name=dorado
#SBATCH --gpus-per-node=A100:1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=124:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output basecallout_%j.out    # save the output into a file
#SBATCH --error basecallerr_%j.err     # save the error output into a file

module purge
module load Dorado/0.5.0

##############DORADO#####################
dorado basecaller sup /nesi/nobackup/uow03920/BlowflyAssemblyData/Promethion_Raw_Data_Blow5/Barcodes_needed_22258/  --recursive --device 'cuda:all' --kit-name SQK-MLK111-96-XL --sample-sheet Barcodes_needed_22258/22258.csv --resume-from /nesi/nobackup/uow03920/BlowflyAssembly/BAM_basecall/22258_sup_calls.bam > /nesi/nobackup/uow03920/BlowflyAssembly/BAM_basecall/22258_sup_calls_rsumd.bam
```

You need to run this twice for each library (if you have two pools), changing the necessary folder paths. The result should be ONE BAM file that contains all the data from the four POD5 files. Next we will demultiplex these so that we have a BAM file for each species. 

# 3. Demultiplexing 

I used the DORADO demux tool (https://github.com/nanoporetech/dorado), running Pool A and Pool B seperately. I tried to run this a couple of times using a slurm script but it failed twice so I just ran it directly on Nesi and it worked fine - took about two hours per library. The output is four BAM files (i.e. one BAM file per sample). 

```
module purge
module load Dorado/0.5.0

cd /nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/01_demux_22258

dorado demux --output-dir nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/01_demux_22258 --no-classify 22258_sup_calls.bam

```

# 4. Convert the BAM files into FASTQ files

Can probably run this directly through NESI without needing to create a slurm script. You need to download samtools - I downloaded from here (http://sourceforge.net/projects/samtools/files/samtools/) and followed these instructions (http://www.sthda.com/english/wiki/install-samtools-on-unix-system). Sometimes needed to re-run this line to get samtools to work (i.e. if nesi had timed out): ```export PATH=$PATH:/nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/samtools-1.20```

We used SAMtools because it is the best for insect genomes. Other tools are highly optimised for humans etc. SAMtools is also the fastest. 

```
cd /nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/01_demux_22258/DEMUX_22258/

samtools fastq SQK-MLK111-96-XL_barcode01.bam > sample1.fastq
samtools fastq SQK-MLK111-96-XL_barcode02.bam > sample2.fastq
samtools fastq SQK-MLK111-96-XL_barcode03.bam > sample3.fastq
samtools fastq SQK-MLK111-96-XL_barcode04.bam > sample4.fastq

cd /nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/02_demux_22276/DEMUX_22276/

samtools fastq SQK-MLK111-96-XL_barcode01.bam > sample1a.fastq
samtools fastq SQK-MLK111-96-XL_barcode02.bam > sample2a.fastq
samtools fastq SQK-MLK111-96-XL_barcode03.bam > sample3a.fastq
samtools fastq SQK-MLK111-96-XL_barcode04.bam > sample4a.fastq
```

# 5. Concatenate the pool A and pool B samples together

Merge the poolA and poolB files together so that there is only one file to work with hereafter

```
cat sample1.fastq  sample1a.fastq > MO_01_cat.fastq
cat sample2.fastq  sample2a.fastq > MO_02_cat.fastq
cat sample3.fastq  sample3a.fastq > MO_03_cat.fastq
cat sample4.fastq  sample4a.fastq > MO_04_cat.fastq

```

# 6. QC of fastQ files using Nanoplot and BBmap

NanoPlot performs quality checks on the FASTQ files.

This is the code for the QC using nanoplot
```
#!/bin/bash -e

#SBATCH --account=uow03920
#SBATCH --job-name=nanoplot
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=124:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output basecallout_%j.out    # save the output into a file
#SBATCH --error basecallerr_%j.err     # save the error output into a file

module purge
module load NanoPlot

#####NANOPLOT#####

for i in 01 02 03 04; do
  NanoPlot --verbose -t 8 --fastq MO_${i}_cat.fastq -o 01_QC/MO_${i}
done
```

BBmap checks the read length statistics (I didn't run this tool because NanoPlot produces these results too but providing the code here just in case I want to use it in the future)

```
for i in 01 02 03 04; do 
readlength.sh in=${i}_cat_fil.fastq out=${i}_histogram.txt
done
```

# 7. Removing duplicate reads 

My FASTQ files had a number of duplicate reads in them so we needed to do an extra step to remove them. This was easy and I just ran it without using a slurm script. P.S at this stage my NeSI account wasn't working properly and I wasn't sure why. Turns out, it was because I had used all of the storage in my nobackup account. SO just make sure that you have space in your nobackup account to do the codes :) 

```
seqkit rmdup MO_01_cat.fastq -n -o MO_01_cat_clean.fastq -D Duplicates.txt -j 16
seqkit rmdup MO_02_cat.fastq -n -o MO_02_cat_clean.fastq -D Duplicates2.txt -j 16
seqkit rmdup MO_03_cat.fastq -n -o MO_03_cat_clean.fastq -D Duplicates3.txt -j 16
seqkit rmdup MO_04_cat.fastq -n -o MO_04_cat_clean.fastq -D Duplicates4.txt -j 16
```
SAMPLE 1 removed 630,142 duplicates;
SAMPLE 2 removed 1,107,312 duplicates;
SAMPLE 3 removed 442,128 duplicates; 
SAMPLE 4 removed 713,723 duplicates.

# 8. Filtering reads using chopper using the 'clean' (i.e. removed duplicates) files

My data had no high molecular weight scores (similar to Meeran's) so we decided to filter my reads based on QUALITY (minimum quality score of 8) and READ LENGTH (minimum length of 500 bases) using Chopper. We decided that --headcrop and --tailcrop were not necessary because Dorado does a good enough job of basecalling. The filtered reads are then saved to new FASTQ files (e.g., MO_${i}_cat_fil.fastq) for each sample.

```
#!/bin/bash -e

#SBATCH --account=uow03920
#SBATCH --job-name=chopper
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=124:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output basecallout_%j.out    # save the output into a file
#SBATCH --error basecallerr_%j.err     # save the error output into a file

module purge
module load chopper

#####CHOPPER#####
for i in 01 02 03 04; do
chopper --threads 8 -q 8 -l 500 < /nesi/nobackup/uow03920/01_Blowfly_Assembly/03_FASTQ/MO_${i}_cat_clean.fastq > /nesi/nobackup/uow03920/01_Blowfly_Assembly/04_Filtered_FASTQ/MO_${i}_cat_clean_fil.fastq ;
done
```

# 9. Repeat quality check again of the new filtered fastq files

Using FASTQC, NanoPlot, BBmap, etc. to check if filtering improved the reads.

# 10. Assemble genomes using FLYE 

I used FLYE to assemble my genomes with a threads score of 16 and an estimated genome size of 700 megabases. We also used a nextflow script to make this work (see below). We used a nextflow script because it is faster than just running through slurm script. This is because nextflow runs each sample parallel (i.e. sample 1 and 2 are being run simultaneously rather than sample 1 being run to completion before sample 2 starts). This assembly using Next Flow took 1 day and 6ish hours. 

NextFlow script (which was named fly_nxtflow.nf)

```
#!/usr/bin/env nextflow

// Define the list of IDs
params.ids = [ '01', '02', '03', '04']

// Create a channel from the list of IDs
Channel.from(params.ids)
    .set { id_ch }

// Define the process for running flye
process runFlye {
    cpus 24
    memory '60 GB'
    publishDir("/nesi/nobackup/uow03920/01_Blowfly_Assembly/04_Filtered_FASTQ/02_FLYE", mode: 'move')
    // Define the input
    input:
    val id from id_ch

    // Define the output
    output:
    file("MO_${id}_flye") into flye_out_ch

    // Define the command to be executed
    script:
    """
    flye --nano-hq /nesi/nobackup/uow03920/01_Blowfly_Assembly/04_Filtered_FASTQ/MO_${id}_cat_clean_fil.fastq \
         --out-dir MO_${id}_flye --threads ${task.cpus} -g 700m
    """
}

// Define what to do with the output (optional)
flye_out_ch.view()
```

Made it run in a slurm script using the following code:

```
#!/bin/bash -e
#SBATCH --account=uow03920
#SBATCH --job-name=nextflow_Flye
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output nextflow_flye_%j.out    # save the output into a file
#SBATCH --error nextflow_flye%j.err     # save the error output into a file

# purge all other modules that may be loaded, and might interfare
module purge

## load tools
module load Flye/2.9.3-gimkl-2022a-Python-3.11.3 
ml Nextflow/22.10.7
### FLYE

nextflow run fly_nxtflow.nf
```



























# Illuminaaaaaaaaa sequenceeeeeees

# Concatenate the two lanes for each sample 

```
cat PI_G_CH_1.fq.gz PI_G_CH_2.fq.gz > PI_G_CH.fq.gz
cat PI_G_CQ_1.fq.gz PI_G_CQ_2.fq.gz > PI_G_CQ.fq.gz
cat PI_G_CS_1.fq.gz PI_G_CS_2.fq.gz > PI_G_CS.fq.gz
cat PI_G_CV_1.fq.gz PI_G_CV_2.fq.gz > PI_G_CV.fq.gz
```

Move all of the concatenated files into a folder together. 

# Check the quality of fastq files using fastqc and multiqc




# Kmer analysis using Illumina reads & Genome quality completeness estimation using Merqury
Follow this tutorial to create kmer a database of high-quality Illumina reads using meryl (https://github.com/marbl/merqury/wiki/1.-Prepare-meryl-dbs). I used a kmer value of 18 rather than the k-mer value I got from running the first line of code included in the tutorial (this was recommended by Meeran and a paper produced by someone from Peter's lab)

``` ml Merqury ``` I include this line of code because I thought I needed to download meryl to do the analysis and spent way too long trying to configure it lol

Do one-step counting on each sample. Each of my samples were in seperate folders, so need to make sure that the directory is the folder with each sample in it

```
meryl k=18k count PI_G_CH_272_DKDL220008012-1A_HHL2KCCX2_L4_1.fq.gz output $genome.meryl
meryl k=18k count PI_G_CQ_50_DKDL220008012-1A_HHL2KCCX2_L4_1.fq.gz output $genome.meryl
meryl k=18k count PI_G_CS_239_DKDL220008012-1A_HHL2KCCX2_L4_1.fq.gz output $genome.meryl
meryl k=18k count PI_G_CV_198_DKDL220008012-1A_HHL2KCCX2_L4_1.fq.gz output $genome.meryl
``` 
