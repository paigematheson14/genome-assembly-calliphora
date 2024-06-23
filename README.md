# Genome assembly - Calliphora sp. 
This repository will detail the workflow I used to create reference genomes for four Calliphora species (C. quadrimaculata, C. hilli, C. stygia, and C. vicina) using PromethION nanopore data and Illumina data. 

# Folder organisation
Creating a genome is a complex process that generates many different files. Thus, it is important to keep folders organised in a way that makes sense. Meeran has a great github page that is dedicated to this and I recommend structuring folders in this manner - https://github.com/meeranhussain/QuickTips

# Sequenced data - PromethION data
My promethion data arrived to me in BLOW5 format (which is funny because blow5 sounds awfully similar to blowfly...). There were two libraries (pool A and pool B) which contained data for four samples (each sample representing one of the four aforementioned species). There were 96 barcodes (in blow5 format) in each library, however, we only really care about the first four barcodes (as these represent our samples and contain the bulk of the data we need). I think that the remaining 92 barcodes just contain junk but I'm not 100% sure on that fact. 

- NB01	GR-AC-CQ-50 (i.e., calliphora quadrimaculata) 
- NB02	GR-AC-CH-259 (i.e., calliphora hilli)
- NB03	GR-AC-CS-235 (i.e., calliphora stygia)
- NB04	GR-AC-CV-196 (i.e., calliphora vicina)

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

The CSV file needs to be in the same folder as the POD5 files. If you have two libraries to basecall, then you need to have two folders and run each library seperately. Here is the slurm script for basecalling: 

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

I used the DORADO demux tool (https://github.com/nanoporetech/dorado), running Pool A and Pool B seperately. Need to create another two slurm scripts to run (one for each library; if you only have one library then you do not need to make two). BTW you only need the '#SBATCH --gpus-per-node=A100:1' line for the basecalling part, so can remove from subsequent slurms :) 

```
#!/bin/bash -e

#SBATCH --account=uow03920
#SBATCH --job-name=dorado
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=8
#SBATCH --time=124:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=paige.matheson14@gmail.com
#SBATCH --output demux_%j.out    # save the output into a file
#SBATCH --error demux_%j.err     # save the error output into a file

module purge
module load Dorado/0.5.0

cd nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/01_demux_22258

dorado demux --output-dir nesi/nobackup/uow03920/01_Blowfly_Assembly/02_basecalling/01_demux_22258 --no-classify 22258_sup_calls.bam

```








