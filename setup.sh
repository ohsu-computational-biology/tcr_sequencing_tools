#!/bin/sh

cd $data/

# Directories to hold raw fastqs and associated files
mkdir fastqs_from_core
mkdir fastqs_from_core/fastqs
mkdir fastqs_from_core/FastQC
mkdir fastqs_from_core/FastQC/html
mkdir fastqs_from_core/FastQC/zip
mkdir fastqs_from_core/FastQC/unzipped
mkdir fastqs_from_core/extras
mkdir fastqs_from_core/md5
mkdir fastqs_from_core/QC_recopy

# Directories to hold mixcr output
mkdir mixcr
mkdir mixcr/alignments
mkdir mixcr/assemblies
mkdir mixcr/despiked_fastqs
mkdir mixcr/export_align
mkdir mixcr/export_clones
mkdir mixcr/indexes
mkdir mixcr/reports
mkdir mixcr/reports/align
mkdir mixcr/reports/assemble

# Directories to hold normalization output
mkdir normalization
mkdir normalization/clones
mkdir normalization/counts
mkdir normalization/normalized_clones
mkdir normalization/QC

# Directories to hold PEAR output
mkdir peared_fastqs
mkdir peared_fastqs/assembled
mkdir peared_fastqs/discarded
mkdir peared_fastqs/unassembled
mkdir peared_fastqs/QC_recopy
mkdir peared_fastqs/QC_recopy/assembled
mkdir peared_fastqs/QC_recopy/discarded
mkdir peared_fastqs/QC_recopy/unassembled

# Directories to hold count spikes output
mkdir spike_counts
mkdir spike_counts/25bp
mkdir spike_counts/25bp/counts
mkdir spike_counts/25bp/qc
mkdir spike_counts/25bp/reads_to_remove
mkdir spike_counts/9bp
mkdir spike_counts/9bp/counts
mkdir spike_counts/9bp/qc
mkdir spike_counts/9bp/reads_to_remove

# QC directory
mkdir QC

# Directory to hold condor's log outputs
mkdir condor_logs
mkdir condor_logs/spike_counts
mkdir condor_logs/spike_counts/25bp
mkdir condor_logs/spike_counts/9bp
mkdir condor_logs/mixcr
mkdir condor_logs/mixcr/align
mkdir condor_logs/mixcr/assemble
mkdir condor_logs/mixcr/despiked
<<<<<<< HEAD
mkdir condor_logs/mixcr/export_align
mkdir condor_logs/mixcr/export_clones
mkdir condor_logs/normalization
<<<<<<< HEAD

# Directory for tools
mkdir tools
mkdir tools/formatted
cp -r $tool/condor_tools/condor_formats $data/tools/
cp -r $tool/condor_tools/submits $data/tools/
=======
=======
mkdir condor_logs/mixcr/export
mkdir condor_logs/normalization
>>>>>>> master
>>>>>>> 221c73a29d8165a7b50df9325fce2e001a6c66fe
