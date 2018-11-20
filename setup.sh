#!/bin/sh

RNA=$1   # Type in 'RNA' to specify that it is an RNA run.

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
mkdir mixcr/clna
mkdir mixcr/vdjca
mkdir mixcr/export_clones
mkdir mixcr/reports

# Different depending on if RNAseq or normal TCRseq
if [ $RNA == "RNA" ]; then
	mkdir mixcr/partial1
	mkdir mixcr/partial2
	mkdir mixcr/extended
	mkdir mixcr/reports/partial1
	mkdir mixcr/reports/partial2
else
	mkdir mixcr/despiked_fastqs
fi



# Directories to hold normalization output
if [ $RNA != "RNA" ]; then
	mkdir normalization
	mkdir normalization/decontam
	mkdir normalization/counts
	mkdir normalization/normalized_clones
	mkdir normalization/QC
fi 

# Directories to hold PEAR output
if [ $RNA != "RNA" ]; then
	mkdir peared_fastqs
	mkdir peared_fastqs/assembled
	mkdir peared_fastqs/discarded
	mkdir peared_fastqs/unassembled
	mkdir peared_fastqs/QC_recopy
	mkdir peared_fastqs/QC_recopy/assembled
	mkdir peared_fastqs/QC_recopy/discarded
	mkdir peared_fastqs/QC_recopy/unassembled
fi

# Directories to hold count spikes output
if [ $RNA != "RNA" ]; then
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
fi

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

# Directory to hold slurm log output
mkdir slurm_logs
mkdir slurm_logs/analysis
mkdir slurm_logs/spike_counts
mkdir slurm_logs/spike_counts/count9
mkdir slurm_logs/spike_counts/count25
mkdir slurm_logs/postprocess
mkdir slurm_logs/postprocess/decontam
mkdir slurm_logs/postprocess/collapse
mkdir slurm_logs/normalize
mkdir slurm_logs/qc
mkdir slurm_logs/remove
mkdir slurm_logs/setup
mkdir slurm_logs/pear
mkdir slurm_logs/mixcr

if [ $RNA == "RNA" ]; then
	mkdir slurm_logs/mixcr/partial1
	mkdir slurm_logs/mixcr/partial2
	mkdir slurm_logs/mixcr/extend
fi

# Directory for tools
mkdir tools
mkdir tools/todo
cp -r $tool/slurm $data/tools/

