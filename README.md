TCRseq_tools
============

This repository contains all of the tools needed to run a TCRseq pipeline from accessing the files on nix (formerly mpssr) to performing QC on output.

GENERAL NOTES
==============

1. Remember to create an interactive job on ExaCloud when running interactive and computationally-intensive jobs.  For example:

     `~% condor_submit -interactive -append 'request_memory = 10 GB' -append 'request_cpus = 4'`  

     `~% srun --mincpus 4 --mem 4G --time 1-00 --pty bash`

The new servers (exahead1 and exahead2) use the slurm job scheduler, not condor, the first command will no longer work.  

DIRECTORY SET UP AND ENV VARS
=============================

1. Log into exacloud via the exahead1 server:

    ```
    ~$ ssh userName@exahead1.ohsu.edu
    ~$ <enter password>
    ```

1. Create a location on ExaCloud for storing the project files, by creating a new directory in the tcrseq project area:

    ```
     ~$ cd /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/
     ~$ mkdir DNAXXXXXLC
     ```

1. Create environment variables for the data and tool directories for easy access to them. You will also need these to be set for the condor\_submit scripts to work correctly.

    a. Open your `~/.bash_profile` (or `~/.bashrc`)

    `~$ emacs .bash\_profile` or `~$ emacs .bashrc`

    b. Enter these lines:

     ```
     data=/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/NEW_DIR/
     export data
     tool=/path/to/this/repo/installation/
     export tool
     ```


    c. Don't forget to source: `~$ . ~/.bashrc`

    d. In the directory you created, create the following folder structure by calling `~$ sh $tool/misc/setup.sh`. If you are running an RNAseq analysis, just add any character string as an argument, in order to get the appropriate directory structre: `~$ sh $tool/misc/setup.sh y`.

```
        .
        ├── condor_logs
        │   ├── mixcr
        │   │   ├── align
        │   │   ├── assemble
        │   │   ├── despiked
        │   │   ├── export_align
        │   │   └── export_clones
        │   ├── normalization
        │   ├── pear
        │   ├── decontaminate
        │   ├── QC
        │   ├── gliph
        │   │   ├── convert
        │   │   └── run
        │   ├── groups
        │   ├── spike_counts
        │   │   ├── 25bp
        │   │   └── 9bp
        │   └── setup
        │   │   ├── md5
        │   │   └── unzip
        ├── fastqs_from_core
        │   ├── extras
        │   ├── FastQC
        │   │   ├── html
        │   │   ├── unzipped
        │   │   └── zip
        │   ├── fastqs
        │   ├── md5
        │   └── QC_recopy
        ├── mixcr
        │   ├── alignments
        │   ├── assemblies
        │   ├── despiked_fastqs
        │   ├── export_align
        │   ├── export_clones
        │   ├── empty_clones
        │   ├── indexes
        │   └── reports
        │       ├── align
        │       └── assemble
        ├── normalization
        │   ├── decontam
        │   ├── counts
        │   ├── normalized_clones
        │   └── QC
        ├── peared_fastqs
        │   ├── assembled
        │   ├── discarded
        │   ├── QC_recopy
        │   │   ├── assembled
        │   │   ├── discarded
        │   │   └── unassembled
        │   └── unassembled
        ├── QC
        ├── spike_counts
        │   ├── 25bp
        │   │   ├── counts
        │   │   ├── qc
        │   │   ├── spikes
        │   │   ├── empty
        │   │   └── reads_to_remove
        │   └── 9bp
        │       ├── counts
        │       ├── qc
        │       └── reads_to_remove
        ├── freqGroups
        │   ├── baseLine
        │   ├── groupData
        │   ├── overLapResults
	│   └── treatSpecificClones
        ├── gliph
        │   ├── byTreat
        │   │   ├── clones
        │   │   └── results
        │   ├── specificTreat
        │   │   ├── clones
        │   │   └── results
	│   └── todo
        └── tools
            ├── condor_formats
            │   ├── 01_format.unzip.R
            │   ├── 02_format.pear.R
            │   ├── 10_format.count.spikes.R
            │   ├── 20_format.remove.spikes.R
            │   ├── 30_format.mixcr.align.R
            │   ├── 40_format.mixcr.assemble.R
            │   ├── 50_format.export.align.R
            │   ├── 51_format.export.pretty.align.R
            │   ├── 60_format.export.clones.R
            │   └── 80_format.normalize.R
            └── submits
                ├── 00_md5.submit
                ├── 01_unzip.submit
                ├── 02_pear.submit
                ├── 10_count.spikes.9bp.submit
                ├── 10_count.spikes.25bp.submit
                ├── 20_remove.spikes.submit
                ├── 30_align.submit
                ├── 40_assemble.submit
                ├── 50_export.align.submit
                ├── 51_export_align_pretty.submit
                ├── 60_export.clones.submit
                ├── 70_decontaminate.submit
                ├── 80_normalize.submit
                ├── 90_runQC.submit
                └── 100_analysis.submit
```

 Unless otherwise stated, instructions below assume that you're at the root of the directory structure stated immediately above.

PREPROCESS
===========

## Copy Files
The Core places the files on their IGL server, nix (formerly mpssr). You will need access to this server in order to copy files to exacloud.

1. Use scp or rsync to copy the fastq files from the Core's IGL server (nix) to the fastqs\_from\_core directory on ExaCloud. For example:  

     ```
     ~% pwd
     /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC
     ~% rsync -r username@mpssr1:/projects/DNAXXXXLC/150*/DNAXXXXLC
     /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/fastqs\_from\_core/fastqs
     ```
1. Use the instructions provided by the Core to find the files on the IGL server.  (You can ssh into the server if necessary – the Core typically provides a temporary password giving you access for ~two weeks.)

1. Copy the fastQC files over to the corresponding directory. The rest of the data can be copied to "extras"

## FastQC
1. Inspect the fastQC files to identify any problematic samples. MultiQC generates an html file that we will have to transfer to our local machine to view. We first have to run MultiQC. To do this, you must either specify the full path to the MultiQC executable, or add `export PATH="/home/exacloud/lustre1/BioCoders/Applications/miniconda3/bin:$PATH"` to your `~/.bashrc` or `~/.bash_profile`. This example assumes you have not added the conda path.

     ```
     ~$  pwd
     /path/to/DNAXXXXLC/fastqs\_from\_core/FastQC/
     ~$ srun --mincpus 1 --mem 16G --time 0-12 --pty bash
     ~$ /home/exacloud/lustre1/BioCoders/Applications/miniconda3/bin/multiqc --interactive --export .
     ```

2. A new file called `multiqc_report.html` will be created in the FastQC directory. You can either copy this to your local, and then view it in a web browser:

     ```
     # THIS IS ON LOCAL COMPUTER # 
     ~$ pwd
     /path/to/temp/dir/
     ~$ scp username@exacloud:/path/to/DNAXXXXLC/fastqs\_from\_core/FastQC/multiqc_report.html ./
     ~$ open multiqc_report.html
     ```

or you can view it directly from exacloud. This option requires that you `ssh` into exacloud with the `-X` flag, and also have XQuartz installed:

```
# LOCAL #
~$ ssh -X username@exacloud
# Remote #
~$ cd /path/to/FastQC/dir/
~$ firefox
# A new firefox browser should automatically open
# Press ctrl + o to search for the multiqc file, and select it
```

## MD5sums
Verify that your file transfer was successful by calculating the MD5sums of the files and comparing them to the MD5sums provided by the core.  

1. The process.md5.R script will compute the MD5sums for the files, which can be checked with those provided by the core:
     ```
     ~% cd ./tools/slurm
     ~% sbatch 00_sbatchmd5.sh
     ```
1. Compare to core values
     ```
     ~$ diff DNA150826LC/calculated.md5.sums.txt DNA150826LC/md5sum.sorted.txt
     ```
1. If any MD5sums differ, retransfer the appropriate files and recheck their MD5sums. Move the three MD5sums files into the md5 directory.

Rename
=======
The naming conventions used by IGL introduce an extra sample tag that is unnesseccary. We want to remove that. 
Optionally, we can remove 1 other field, if there happens to be a treatment designation or something that is also included.

```
~$ Rscript $tool/10_preProcess/01_rename_fastqs.R /path/to/fastqs/ /path/to/QC/
```

UnZip
=====
The files are zipped to allow fast transfer from nix, but need to be unzipped for the remainder of the pipeline.  
We will now begin using the slurm job scheduler to run our jobs. For each step, there is generally a .sh or .R file that runs the program,
and an sbatch script that will submit the program to the job scheduler. A few modifications to the sbatch scripts will need to be made.

```
~$ cd $data/tools/slurm
~$ vi 02_sbatchUnzip.sh
< change --array argument to be 0-X, where X is the largest tens place in your files. 130 files would be 0-13 >
~$ sbatch 01_sbatchUnzip.sh
~$ mv unzip\* $data/slurm_logs/setup/ 
```

After each step, move the .out and .err files to the appropriate log directory.
For subsequent steps, 'todo' files will be created, which is described below. If you are running an RNAseq pipeline, you can skip the PEAR and PROCESS SPIKES steps.

PEAR
=====
We receive paired-end sequencing files from IGL. Before we can use these files in later steps, we need to merge the forward and reverse reads. 
We do this using PEAR.

```
~$ cd $data/fastqs_from_core/fastqs
~$ ls -v *_R1_* > $data/tools/todo/pear.txt
~$ cd $data/tools/slurm/
~$ vi 03_sbatchPear.sh
< change --array argument to be 1-X, where X is the number of files. 130 files would be 1-130 >
~$ sbatch 03_sbatchPear.sh
```

Sometimes the job scheduler is pretty busy and running an entire batch at once causes some of the jobs to fail.  
Instead of running them all at the same time, you can utilize the `runSbatch.sh` script. This script will submit jobs in batches of 10,
waiting for one set of 10 to finish before submitting another. If there were 130 samples in a batch, it would look like:

```
~$ sbatch runSbatch.sh 03_sbatchPear.sh 130
```

The PEAR program produces fastq files of the assembled reads, but it also produces files containing all discarded reads, 
as well as all unassembled reads. Everything is output the the `$data/peared_fastqs/assembled/` directory, 
so the other files will need to be moved to their corresponding directories.

PROCESS SPIKES
===============
## Count Spikes
From now on, we will only be using the assembled fastq files in the analysis. Our reads have synthetic spike-ins along with the TCR DNA. 
We use a 9-bp barcode to remove all spikes, and a 25-bp barcode to quantify the levels of each spike.

```
~$ cd $data/peared_fastqs/assembled/
~$ ls -v * > $data/tools/todo/count9.txt
~$ cd $data/tools/slurm
~$ vi 10_sbatchCountSpikes9.sh
< change --array argument to be 1-X >
~$ sbatch 10_sbatchCountSpikes9.sh
```
Run this for both 9bp and 25bp.  

## Remove spikes
We need to remove all of the spike-in sequences from our fastq files so that we only submit true TCR sequences to MiXCR. 
This uses the same base input as count spikes, so you can do: `cp $data/tools/todo/count9.txt $data/tools/todo/despike.txt`. Then just
open `20_sbatchRemoveSpikes.sh` and change the `--array` variable, as previously described, and submit. 

IDENTIFY CLONOTYPES
===================
Our reads are finally ready to be processed through the MiXCR analysis pipeline. This process differs depending on if the data are DNA or RNA. 

## DNA  
Now with MiXCR version 3.0+, there is only one tool, `analyze`, which will run all of the old `align`, `assemble`, and `export` modules. Either create the `todo/analyze.txt` file using the `mixcr/despiked_fastqs` directory, or use the non-todo method. See below for the documentation for MiXCR version 2.x.  

**Old version**  
We use three tools from the MiXCR suite: `align`, `assemble`, and `export`. Generally use the output of the previous command to form the "todo" file of the next command. 
So the files in `despiked` go into `todo/align.txt`, the files in `align` will go into `todo/assemble.txt`, and the files in `assemble` will go into `todo/exportClones.txt`.
The tool `50_sbatchMixcrExportAlign.sh` is not used in the standard pipeline and only needs to be run when troubleshooting. For each tool, create the todo file, edit the array argument, submit the sbatch,
move the log files, and continue to the next tool. 

```
~$ 30_sbatchMixcrAlign.sh <--------> todo/align.txt
~$ 40_sbatchMixcrAssemble.sh <-----> todo/assemble.txt
~$ 60_sbatchMixcrExportClones.sh <-> todo/exportClones.txt
```

You must complete `align` before you can run `assemble`.  
Sometimes a few of the clonotype files do not properly form. To check that `exportClones` has run successfully, 
run: `sh $tool/misc/check_clones.sh /path/to/export_clones`. If any files are output, 
go to your submission file and re-run only those files by commenting out the rest.

## RNA
**Note** this has not been updated for MiXCR 3.0+!!
Due to the nature of randomly digested RNAseq data, a few extra steps need to be taken to reliably identify clones. The RNAseq-specific submission files are in
`$data/tools/slurm/rnaseq`. The same principle applies, however, of creating a todo file from the previous step's output, editing the array argument for the appropriate
number of samples, submitting, and moving the log files. The steps for this are:

```
~$ 31_sbatchMixcrAlignRNAseq.sh <------> todo/align.txt
~$ 41_sbatchMixcrAssemblePartial1.sh <-> todo/partial1.txt
~$ 42_sbatchMixcrAssemblePartial2.sh <-> todo/partial2.txt
~$ 51_sbatchMixcrExtendAlign.sh <------> todo/extend.txt
~$ 55_sbatchMixcrAssemble.sh <---------> todo/assemble.txt
~$ 60_sbatchMixcrExportClones.sh <-----> todo/exportClones.txt
```

DECONTAMINATE
===============
After one of our control sequencing runs, some monoclonal sequences began contaminating our reads. We want to remove them from our files before we do anything else. 
This script works on the entire directory, so there is no todo file. Additionally, downstream steps depend on this format, so be sure to run this even if you don't expect
any contamination.  

```
~$ cd $data/tools/slurm
~$ sbatch 40_sbatchDecontaminate.sh
```

NORMALIZE
==========  
Skip this step if you have RNAseq data.  
We have to do some more pre-processing before we are ready to normalize. 
1. Copy the 25-bp spike files to the normalization directory.
1. Calculate the scaling factor
     ```
     ~$ cd /path/to/normalization/dir/
     ~$ Rscript $tool/40_postProcess/calculate.scaling.factor.R /path/to/normalization/counts/ $tool/00_reference/text_barcodesvj.txt 
     ```
1. Normalize the clones
     ```
     ~$ cd $data/tools/slurm
     ~$ vi 50_sbatchNormalize.sh
       < change --array argument to be 1-X >
     ~$ sbatch 50_sbatchNormalize.sh
     ```

POSTPROCESS
=============
Now we need to run all of the QC scripts. Each script can be run individually, 
but they have been placed inside of a bash script that will run them all simultaneously.
```
~$ cd $data/tools/slurm
~$ sbatch 60_sbatchQC.sh
```
ANALYZE
========
There are a few different standard analyses that can be run. 
Some of these require a metadata file that contains treatment information.
The one that is always run is `70_sbatchDiversityAnalysis.sh` 
which generates a text file containing standard diversity statistics such as clonality, shannon entropy, etc. 
There are also variations of this script that allow you to subset for clones 
of certain frequency groups only or certain numbers of top clones.  

There is also the clonotype homeostasis analysis, 
`80_sbatchClonalDivisionSummary.sh` which creates a few excel workbooks and some homeostasis plots. 
Finally you can group clones and subset them based on those groupings.  

Run any of these scripts just like you have been for the other sbatch scripts.

Summary Output
===============
After running all of the QC scripts, and the analysis, combine the outputs into an excel workbook. 
Run the combineQC.R script in the QC tool directory:

```
~$ Rscript $tool/QC/combineQC.R /path/to/QC/dir /path/to/QC/dir
```

You should have the following sheets, comprised of the corresponding QC file:
```
 PEAR:   pear\_summary\_log.txt
 counts:  count.spikes.QC.summary.txt
 spike\_removal:  remove.spikes.QC.result.txt
 mixcr\_alignment:  mixcr.alignment.QC.summary.txt
 mixcr\_assembly:  mixcr.assemble.QC.summary.txt
 decontamination: DNAXXXXLC_contaminationQC.txt
 normalization: aggregate\_normalization\_factor\_QC.txt
 nonStandard_VHits: nonStandard_VHits.txt
 analysis: uniques.shannon.clonality.txt
```


