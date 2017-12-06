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
mkdir mixcr/empty_clones
mkdir mixcr/indexes
mkdir mixcr/reports
mkdir mixcr/reports/align
mkdir mixcr/reports/assemble

# Directories to hold normalization output
mkdir normalization
mkdir normalization/decontam
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
mkdir spike_counts/9bp/spikes
mkdir spike_counts/9bp/empty

# Directories to hold freqGroups output
mkdir freqGroups
mkdir freqGroups/baseLine
mkdir freqGroups/groupData
mkdir freqGroups/overLapResults
mkdir freqGroups/treatSpecificClones

# Directories to hold gliph output
mkdir gliph
mkdir gliph/byTreat
mkdir gliph/byTreat/clones
mkdir gliph/byTreat/results
mkdir gliph/specificTreat
mkdir gliph/specificTreat/clones
mkdir gliph/specificTreat/results

# QC directory
mkdir QC
mkdir QC/meta
mkdir QC/std
mkdir QC/homeo

# Directory to hold condor's log outputs
mkdir condor_logs
mkdir condor_logs/spike_counts
mkdir condor_logs/spike_counts/25bp
mkdir condor_logs/spike_counts/9bp
mkdir condor_logs/mixcr
mkdir condor_logs/mixcr/align
mkdir condor_logs/mixcr/assemble
mkdir condor_logs/mixcr/despiked
mkdir condor_logs/mixcr/export_align
mkdir condor_logs/mixcr/export_clones
mkdir condor_logs/normalization
mkdir condor_logs/pear
mkdir condor_logs/decontaminate
mkdir condor_logs/QC
mkdir condor_logs/setup
mkdir condor_logs/setup/unzip
mkdir condor_logs/setup/md5
mkdir condor_logs/aggregate
mkdir condor_logs/gliph
mkdir condor_logs/gliph/run
mkdir condor_logs/gliph/convert
mkdir condor_logs/groups

# Directory to hold slurm log output
mkdir slurm_logs
mkdir slurm_logs/mixcr
mkdir slurm_logs/mixcr/align
mkdir slurm_logs/analysis
mkdir slurm_logs/mixcr/assemble
mkdir slurm_logs/spike_counts
mkdir slurm_logs/spike_counts/count9
mkdir slurm_logs/spike_counts/count25
mkdir slurm_logs/postprocess
mkdir slurm_logs/postprocess/decontam
mkdir slurm_logs/postprocess/collapse
mkdir slurm_logs/mixcr/exportAlign
mkdir slurm_logs/mixcr/exportClones
mkdir slurm_logs/normalize
mkdir slurm_logs/qc
mkdir slurm_logs/remove

# Directory for tools
mkdir tools
mkdir tools/todo
cp -r $tool/condor_tools/condor_formats $data/tools/
cp -r $tool/condor_tools/submits $data/tools/
cp -r $tool/slurm $data/tools/

