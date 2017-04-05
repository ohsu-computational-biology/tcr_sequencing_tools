TCRseq_tools
============

This repository contains all of the tools needed to run a TCRseq pipeline from accessing the files on mpssr to performing QC on output.

GENERAL NOTES
==============

1. Remember to create an interactive job on ExaCloud when running interactive and computationally-intensive jobs.  For example:

     `~% condor_submit -interactive -append 'request_memory = 10 GB' -append 'request_cpus = 4'`

DIRECTORY SET UP AND ENV VARS
=============================

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

    d. In the directory you created, create the following folder structure by calling `~$ sh $tool/setup.sh`

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
The Core places the files on their IGL server, mpssr. You will need access to this server in order to copy files to exacloud.

1. Use scp or rsync to copy the fastq files from the Core's IGL server (mpssr) to the fastqs\_from\_core directory on ExaCloud.  For example:  

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
     ~$ condor_submit -interactive -append "request_memory = 32 GB"
     ~$ /home/exacloud/lustre1/BioCoders/Applications/miniconda3/bin/multiqc .
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
     ~% cd ./tools/submits
     ~% condor_submit 00_md5.submit
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
~$ Rscript $tool/misc/rename_fastqs.R /path/to/fastqs/ /path/to/QC/
```

UnZip
=====
The files are zipped to allow fast transfer from mpssr1, but need to be unzipped for the remainder of the pipeline. 
We will now begin using the condor job scheduler to run our jobs. This basic format will be followed for most of the remaining steps. 
All the format scripts point to 1-3 input locations and the final argument should almost always be `../submits` for the output of the final submission script. 
You can use `head -15 format.file.R` to view the top of the file and see example locations for the inputs.

```
~$ cd /path/to/tools/condor_formats/
~$ Rscript 01_format.unzip.R /path/to/fastq/dir ../submits
~$ cd ../submits
~$ condor_submit 01_unzip.submit
```
For subsequent steps, open each format file to see what arguments it needs. Then follow this basic procedure.

PEAR
=====
We receive paired-end sequencing files from IGL. Before we can use these files in later steps, we need to merge the forward and reverse reads. 
We do this using PEAR.
```
~$ cd /path/to/tools/condor_formats/
~$ Rscript 02_format.pear.R /path/to/fastq/dir ../submits
etc.
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
~$ cd /path/to/tools/codor_formats/
~$ Rscript 10_format.count.spikes.R /path/to/assembled/fastqs/ [#bp] ../submits
etc.
```
Don't forget the [#bp] argument. It should either be 9 or 25. This needs to be run once for each. 

## Remove spikes
We need to remove all of the spike-in sequences from our fastq files so that we only submit true TCR sequences to MiXCR. 

```
~$ cd /path/to/tools/condor_formats/
~$ Rscript 20_format.remove.spikes.R /path/to/assembled/fastqs/ /path/to/9bp/reads.to.remove/ ../submits
etc.
```

IDENTIFY CLONOTYPES
===================
Our reads are finally ready to be processed through the MiXCR analysis pipeline. 
We use three tools from the MiXCR suite: `align`, `assemble`, and `export`.

```
~$ cd /path/to/tools/condor_formats/
~$ Rscript 30_format.mixcr.align.R /path/to/mixcr/despiked/ ../submits
etc.
```
You must complete `align` before you can run `assemble`. Once both of these are completed, you can run `exportAlign` and `exportClones` simultaneously. 
Sometimes a few of the clonotype files do not properly form. To check that `exportClones` has run successfully, 
run: `sh $tool/misc/check_clones.sh /path/to/export_clones`. If any files are output, 
go to your submission file and re-run only those files by commenting out the rest.

DECONTAMINATE
===============
After one of our control sequencing runs, some monoclonal sequences began contaminating our reads. We want to remove them from our files before we do anything else. 
This script works on the entire directory, and so we don't need to format a submission script:
```
~$ cd /path/to/tools/submits
~$ condor_submit 70_decontaminate.submit
```

NORMALIZE
==========
We have to do some more pre-processing before we are ready to normalize. 
1. Copy the 25-bp spike files to the normalization directory.
1. Calculate the scaling factor
     ```
     ~$ cd /path/to/normalization/dir/
     ~$ Rscript $tool/normalize/calculate.scaling.factor.R /path/to/normalization/counts/ $tool/reference/text_barcodesvj.txt 
     ```
1. Normalize the clones
     ```
     ~$ cd /path/to/tools/submits
     ~$ Rscript 80_format.normalize.R /path/to/normalization/decontam /path/to/normalization/counts/ ../submits
     etc.
     ```

POSTPROCESS
=============
Now we need to run all of the QC scripts. Each script can be run individually, 
but they have been placed inside of a bash script that will run them all simultaneously.
```
~$ cd /path/to/tools/submits
~$ condor_submit 90_runQC.submit
```
ANALYZE
========
Calculate clonotype diversity statistics (clonality, shannon entropy, etc.) using the analysis.R script. 
There is a submission file in the `condor_tools/submissions` directory.

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


