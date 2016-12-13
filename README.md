# TCRseq_tools

All of the tools needed to run a TCRseq pipeline from accessing the files on mpssr to performing QC on output.



**GENERAL NOTES**

1. Remember to create an interactive job on the ExaCloud when running interactive and computationally-intensive jobs.  For example:

     `~% condor_submit -interactive -append 'request_memory = 10 GB' -append 'request_cpus = 4'`

2. The Core places the files on their IGL server, mpssr. You will need access to this server in order to copy files to exacloud.


**DIRECTORY SET UP AND ENV VARS**

1. Create a location on ExaCloud for storing the fastq files, by creating a new directory here:

    `/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/`

2. Create environment variables for the data and tcr\_sequencing\_tools for easy access to them. You will also need these to be set for the condor\_submit scripts to work correctly.

    a. Open your bash profile (or bashrc)

    `~$ emacs .bash\_profile` or `~$ emacs .bashrc`

    b. Enter these lines:

        data=/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/NEW_DIR/
        export data
        tool=/home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/
        export tool

    c. In the directory you created, create the following folder structure by calling:

    `~$ sh $tool/setup.sh`
    
        .
        ├── condor_logs
        │   ├── mixcr
        │   │   ├── align
        │   │   ├── assemble
        │   │   ├── despiked
        │   │   ├── export_align
        │   │   └── export_clones
        │   ├── normalization
        │   ├── pear
        │   ├── spike_counts
        │   │   ├── 25bp
        │   │   └── 9bp
        │   └── unzip
        ├── fastqs_from_core
        │   ├── extras
        │   ├── FastQC
        │   │   ├── html
        │   │   ├── unzipped
        │   │   └── zip
        │   ├── fastqs
        │   ├── md5
        │   └── QC_recopy
        ├── mixcr
        │   ├── alignments
        │   ├── assemblies
        │   ├── despiked_fastqs
        │   ├── export_align
        │   ├── export_clones
        │   ├── indexes
        │   └── reports
        │       ├── align
        │       └── assemble
        ├── normalization
        │   ├── clones
        │   ├── counts
        │   ├── normalized_clones
        │   └── QC
        ├── peared_fastqs
        │   ├── assembled
        │   ├── discarded
        │   ├── QC_recopy
        │   │   ├── assembled
        │   │   ├── discarded
        │   │   └── unassembled
        │   └── unassembled
        ├── QC
        ├── spike_counts
        │   ├── 25bp
        │   │   ├── counts
        │   │   ├── qc
        │   │   └── reads_to_remove
        │   └── 9bp
        │       ├── counts
        │       ├── qc
        │       └── reads_to_remove
        └── tools
            ├── condor_formats
            │   ├── 00.1_format.unzip.R
            │   ├── 00_format.pear.R
            │   ├── 01_format.count.spikes.R
            │   ├── 02_format.remove.spikes.R
            │   ├── 03_format.mixcr.align.R
            │   ├── 04_format.mixcr.assemble.R
            │   ├── 05_format.export.align.R
            │   ├── 06_format.export.clones.R
            │   └── 07_format.normalize.R
            ├── formatted
            └── submits
                ├── 00.1_unzip.submit
                ├── 00_pear.submit
                ├── 01_count.spikes.submit
                ├── 02_remove.spikes.submit
                ├── 03_align.submit
                ├── 04_assemble.submit
                ├── 05_export.align.submit
                ├── 06_export.clones.submit
                └── 07_normalize.submit


 
Unless otherwise stated, instructions below assume that you&#39;re at the root of the directory structure stated immediately above.

**PREPROCESS**

Use scp to copy the fastq files from the Core&#39;s IGL server (mpssr) to the &quot;fastqs\_from\_core&quot; directory on ExaCloud.  For example:

        ~% pwd

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA150826LC

        ~% nohup scp -r leyshock@igl-fs.ohsu.edu:/projects/DNA150826LC/150\*

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA150826LC/fastqs\_from\_core/fastqs

Use the instructions provided by the Core to find the files on the IGL server.  (You can ssh into the server if necessary – the Core typically provides a temporary password giving you access for ~two weeks.)

Copy the fastQC files over as well.

Inspect the fastQC files to identify any problematic samples.

Scp the zip fastQC files to a temporary directory on your local machine and run MultiQC

        ~$ pwd

                /Users/local\_machine/

~$ scp username@exacloud.oshu.edu:/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/fastqs\_from\_core/FastQC/zip/\* ~/Desktop/DNAXXXXLC/FastQC

        ~$ multiqc .

Take the multiqc\_report\_data files and scp them into the QC directory on exacloud.

Verify that your file transfer was successful by calculating the MD5sums of the files and comparing them to the MD5sums provided by the core.  Use the process.md5.R script if useful; for example, compute the MD5sums for the files:

 ~% Rscript /home/exacloud/lustre1/lCompBio/genomic\_resources/tcrseq/md5\_processing/process.md5.R &quot;DNA150826LC/&quot;

then compare them to the Core&#39;s values:

        ~$ diff DNA150826LC/calculated.md5.sums.txt DNA150826LC/md5sum.sorted.txt

If any MD5sums differ, retransfer the appropriate files and recheck their MD5sums. Move the three MD5sums files into the md5 directory.

Unzip the fastq files:

        ~$ pwd

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA150826LC/fastqs\_from\_core/fastqs

        ~$ gunzip \*

Rename the files, eliminating the 1-5 character treatment tag that follows DNAXXXLC\_ and placing it in a vector for use in vdjtools later.

        ~$ pwd

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA160107LC/fastqs\_from\_core/fastqs/

        ~$ Rscript /home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/rename\_fastqs.R ./

Use PEAR to merge the forward and reverse fastq files.  Place results in the appropriate directory:

 peared\_fastqs

   assembled

   discarded

   unassembled

For example:

~$ perl /home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/PEAR/scripts/run\_pear\_exacloud.pl

/home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA150624LC/fastqs\_from\_core/fastqs/\*.fastq

        -o /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA150624LC /peared\_fastqs

        -p 4

Copy pear\_summary\_log.txt to the QC directory.

**PROCESS SPIKES**

Moving forward, we only use the assembled fastq files in:

 peared\_fastqs

  assembled

Count the number of 9bp and 25bp spikes in the fastqs.  To do so, create the appropriate condor.submit file.  First copy the condor\_tools directory into the data directory:

        ~$ pwd                /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/

        ~$ cp –r $tool/tcr\_sequencing\_tools/condor\_tools/\*

Add a directory called &quot;formatted&quot; and another called &quot;submits&quot;.

To generate the repeating chunks of code in the condor.submit file, use the format.count.spikes.condor.R script, located in the condor\_tools directory:

        ~$ pwd

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNAXXXXLC/condor\_tools/

        ~$ Rscript format.sount.spikes.condor.R /fastq/dir/path bp direction

When you run the format.for.condor() function defined in the script, be sure to specify the number of base pairs (either &quot;9bp&quot; or &quot;25bp&quot;, typically) and the direction of search (either &quot;fwd&quot; or &quot;rev&quot;, typically). Open the script to see an example of the correct dir path. Move the output file to the formatted directory.

Copy the submit template and rename it to a .submit file

 ~$ cp ./submit\_templates/count.spikes.template ./submits/25bp.count.spikes.submit

Paste the output of the formatting script into the .submit file. You shouldn&#39;t need to edit any values in the final submit script, but run a test case with the 1

# st
 file to check.

Once the condor.submit file is prepared, submit the jobs to Condor:

 ~$ condor\_submit /home/users/leyshock/condor\_submits/25bp\_forward\_parallel.submit

Place the outputs in the appropriate directory:

 spike\_counts

  25bp

   counts

   qc

   reads\_to\_remove

and place the QC outputs in the QC directory.

Repeat the process above, this time looking for 9bp spikes.  Place the outputs in the appropriate directory:

 spike\_counts

  9bp

   counts

   qc

   reads\_to\_remove

and place the QC outputs in the QC directory.

Remove spiked records from the fastq files. First, create an empty reverse file in the reads\_to\_remove directory.

        ~$ pwd

        /home/exacloud/lustre1/CompBio/data/tcrseq/dhaarini/DNA151124LC/spike\_counts/9bp/reads\_to\_remove

        ~$ touch reverse.txt

 Then, create the appropriate condor.submit file.  Use the files here as templates:

/home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/tcr\_sequencing\_tools/condor\_tools/submit\_templates

To generate the repeating chunks of code in the condor.submit file, use the format.remove.spikes.condor.R script, located here:

 ~$ pwd

/home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/tcr\_sequencing\_tools/condor\_tools

~$ Rscript format.remove.spikes.condor.R /path/to/peared/fastqs /path/to/reads/to/remove

Note that you need to pass multiple collections of files to format.for.condor(), including the fastq files, and one or more reads.to.remove.txt files. Open the script to see examples of correct paths.

Rename the template as .submit instead of .template. Paste the output of the formatting script into the .submit file. Modify the appropriate variables in the &quot;header&quot; of the condor.submit file, to indicate the location of your files, location you&#39;d like Condor&#39;s logger to write to, etc.

Once the condor.submit file is prepared, submit the jobs to Condor:

 ~$ condor\_submit /home/users/leyshock/condor\_submits/ XYZ.submit



Place the despiked fastqs here:

mixcr

  despiked\_fastqs

and place the QC outputs in the QC directory.

**IDENTIFY CLONOTYPES**

Run the reads with spikes removed through MiXCR.  We use three of MiXCR&#39;s functions; in order:  align(), assemble(), and exportClones().  For each step, you&#39;ll need to create the appropriate condor.submit file.  Use the files here as templates:

/home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/tcr\_sequencing\_tools/condor\_tools/submit\_templates/

To generate the repeating chunks of code in the condor.submit file, use the appropriate script, located here:

 ~$ pwd

/home/exacloud/lustre1/CompBio/genomic\_resources/tcrseq/tcr\_sequencing\_tools/condor\_tools/

~$ Rscript format.mixr.[function].R options

You&#39;ll see three scripts there, one for each MiXCR function:

format.mixcr.align.R

format.mixcr.assemble.R

format.mixcr.export.R

Open each script to see examples of directory paths to use as command line options. Rename the .template file into a .submit file. Paste the output of the formatting script into the .submit file.  Modify the appropriate variables in the &quot;header&quot; of the condor.submit file, to indicate the location of your files, location you&#39;d like Condor&#39;s logger to write to, etc.

Once the condor.submit file is prepared, submit the jobs to Condor:

 ~$ condor\_submit /home/users/leyshock/condor\_submits/ XYZ.submit



**NORMALIZE**

Prepare the data for normalization, by arranging your directories.  Copy the clonotypes you exported from MiXCR here:

 normalization

  clones

and copy the 25bp spike count files (e.g. &quot;S27.assembled.spike.counts.25bp.txt&quot;) here:

 normalization

  counts

Calculate the normalization scaling factor by running the calculate.scaling.factor.R R script (in tcr\_sequencing\_tools/normalization/).  Pass as an input to the function the &quot;counts&quot; directory:

 normalization

  counts

Normalize the clones with the normalize.clonotype.counts.condor.R script.  Create the appropriate condor.submit file.  Use the files here as an example:

 /mnt/lustre1/CompBio/genomic\_resources/tcrseq/condor\_examples

To generate the repeating chunks of code in the condor.submit file, use the &quot;format.normalize.condor.R&quot; script, located here:

 /mnt/lustre1/CompBio/genomic\_resources/tcrseq/tcr\_sequencing\_tools/condor\_tools

Note that you&#39;ll need to pass the format.for.condor() function multiple arguments.  Paste the output of the formatting script into the condor.submit file.  Modify the appropriate variables in the &quot;header&quot; of the condor.submit file, to indicate the location of your files, location you&#39;d like Condor&#39;s logger to write to, etc.

Once the condor.submit file is prepared, submit the jobs to Condor:

 ~$ condor\_submit /home/users/leyshock/condor\_submits/XYZ.submit

Place the normalized clonotype files in a new directory:

 normalization

  normalized\_clones

**POSTPROCESS**

Collect the various QC outputs and place in the &quot;QC&quot; directory.  Typically the QC outputs include &quot;QC&quot; or &quot;qc&quot; in the file name.

Create a temporary subdirectory and place the normalization QC outputs (one for each sample) there.  Run the &quot;aggregate.normalization.QC.R&quot; script, passing the temp directory as a parameter to the script.  Place the output file in the root QC directory.

Clear out the temp directory, then place all of the count.spike files there; they&#39;ll have a structure similar to &quot;sample\_id.assembled.qc.txt&quot;.  Run the &quot;count.spikes.QC.R&quot; script, passing the directory as an input to the script.

Continue to aggregate the QC results for both the MiXCR alignment and assembly steps.  This will generate, respectively, files &quot;mixcr.alignment.QC.summary.txt&quot; and &quot;mixcr.assemble.QC.summary.txt&quot;.

Import all of the QC files into an Excel workbook, and name the worksheets appropriately.  Your final QC workbook should include these tabs (with respective inputs):

 PEAR:   pear\_summary\_log.txt

 counts:  count.spikes.QC.summary.txt

 spike\_removal:  remove.spikes.QC.result.txt

 mixcr\_alignment:  mixcr.alignment.QC.summary.txt

 mixcr\_assembly:  mixcr.assemble.QC.summary.txt

 mixcr\_export: (not created yet)

 normalization: aggregate\_normalization\_factor\_QC.txt



**ANALYZE**

Calculate clonotype entropies by running the script uniques.shannon.clonality.R. This also calculates the clonality, which is a different measure of diversity, as well as the total number of unique clones. If the batch requires it, also run compare.replicates.R.

**FINAL STEPS**

Aggregate the QC workbook with the analysis output into 1 final workbook. In addition, add the metadata.txt file from vdjtools (or DM&#39;s sample identification file). And upload to the appropriate directory in Box.