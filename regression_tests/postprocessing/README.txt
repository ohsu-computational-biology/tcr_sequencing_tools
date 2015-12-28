
The files in this directory can be used to test:

1.  frameshift.and.stop.remover.R

To test, start an R session, and load the file above.  Depending on whether or not
    remove.clone.table is set to TRUE or FALSE (default is FALSE), one of two 
    files will be output.

Results will be output to the same directory as the input file, but the file will
    have "_postprocessed.txt" appended to it.  

Diff the results output with the appropriate file in the expected_outputs/ 
    subdirectory:

        - ".no.remove" (records not removed)
        - ".removed" (records removed)


