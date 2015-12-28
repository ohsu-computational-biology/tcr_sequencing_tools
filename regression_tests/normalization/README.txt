
The files in this directory can be used to test:

1.  normalize.clonotype.counts.condor.R

Run this command:

	~$ Rscript /Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/normalize/normalize.clonotype.counts.condor.R /Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/regression_tests/normalization/inputs/clones/S1_clones_exported.txt /Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/regression_tests/normalization/inputs/counts/S1.spike.counts.txt /Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/normalization/ /Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/regression_tests/normalization/inputs/scaling_factor.txt 

the results will be output to:

	/Users/leyshock/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/regression_tests/normalization/

Diff the results output with this file:

	/Users/leyshock/Desktop/TCRseq/tools/tcr_sequencing_tools/regression_tests/normalization/expected_outputs/S1.exported.clones.normalized.txt

Mutatis mutandis, perform the same operations to test S2.

