"""
fastqcount.py: counts occurrences of specified spikes in a fastq file
Matthew Jagielski - jagielsk@ohsu.edu
"""

# import everything
import argparse # for parsing arguments
import sys
from os import path, walk # some filesystem functions
from collections import defaultdict
from random import randint
#import ctypes

def parse_reads(path,sharedstrs,spikedata,spikeli,maxerrs,removal):
    """
    parse_reads: reads the file, then calls process on each read, managing the list of spikes
    Arguments:
        path: the path to the fastq file
        sharedstrs, spikedata: see the process_spikes docstring for details
        spikeli: the list of spikes, see parse_spikes docstring for details
        maxerrs: the maximum number of allowed errors in a spike transcription
        removal: True if spiked reads should be removed
    Returns:
        spikedict: a dictionary with keys as spikes (SPIKE_ID, SPIKE) and values as counts of the spike
    """

    reverse = path.endswith('R2_001')
    spikedict={}
    for spike in spikeli[reverse][1:]:
        spikedict[spike]=0
    spikedict[spikeli[0][0]] = "COUNT"
    spikedict[('0','')]=0 # initialize
    count =0
    read_id_file = open(path+'.reads.to.remove.txt', 'w')
    with open(path+'.fastq') as genefile:
        if removal:
            geneoutfile = open(path+'rm.fastq','w')
            currecord=[]
            toremove=False
            for i,record in enumerate(genefile):
                currecord.append(record)
#                if i%10000==0: # progress marker
#                    print i
                if i%4==1:
                    curspike = process(record.strip(),sharedstrs[reverse],spikedata,spikeli[reverse],maxerrs)
                    if curspike!=('0',''):
                        toremove=True
                    spikedict[curspike]+=1
                if i%4==3:
                    if not toremove:
                        for line in currecord:
                            geneoutfile.write(line)
                    currecord=[]
                    toremove=False
            geneoutfile.close()
        else:
            i=0
	    current_read_id = genefile.readline()
            while current_read_id:
                i+=1
#                if i%10000==0:
#                    print i
		search_result = process(genefile.readline(),sharedstrs[reverse],spikedata,spikeli[reverse],maxerrs)
		#	Check for mismatch between returned value and expected
		#	This "if" comparison is ugly but it works.  TODO:  do the comparison against the tuple
		ssr = "".join(search_result)
		if ssr != "0":
#			print current_read_id
			read_id_file.write(current_read_id)
		spikedict[search_result]+=1
                genefile.readline()
                genefile.readline()
		current_read_id = genefile.readline()
    return spikedict
            
def process(read,sharedstrs,spikedata,spikeli,maxerrs):
    """
    process: finds spikes in a read
    Arguments:
        read: the read to be searched for spikes
        sharedstrs, spikedata: see process_spikes docstring for details
        spikeli: the list of spikes, see parse_spikes docstring for details
        maxerrs: the maximum number of allowed errors in the sequence
    Returns:
        spike: the spike found in the read, or a blank spike if no spike exists
    """
    spikelen=34 # spikes have 34 characters - change this if that is wrong
    m = len(sharedstrs[0]) # should be 9
    for i in range(len(read)-spikelen+1):
        requiredchars=True # check that this starting position has all required characters
        for charno in range(m):
            if read[i+charno]!=sharedstrs[0][charno]:
                requiredchars=False
                break
            if read[i+spikelen-charno-1]!=sharedstrs[1][-1-charno]:
                requiredchars=False
                break
        if requiredchars: # check to see if there is a spike iff all required characters are present
            #print read[i:i+spikelen]
            for spike in spikeli:
                if read[i:i+spikelen]==spike[1]:
                    return spike
    return ('0','') # the default value

def parse_spikes(path):
    """
    parse_spikes: reads spike configuration file into a list
    Arguments: path - path to the spike configuration file
    Returns: spikeli - the list of spikes in the following format: (SPIKE_ID, SPIKE)
    """
    spikeli = [[],[]]
    with open(path) as spikefile:
        for line in spikefile:
            vals=line.split()
            spikeli[0].append(tuple(vals)) # add all spikes to the list as (SPIKE_ID, SPIKE)
            vals[1] = conjugate(vals[1])
            spikeli[1].append(tuple(vals))
    return spikeli
    
def process_spikes(spikeli):
    """
    process_spikes: converts list of spikes into more detailed data about all spikes
    Arguments: spikeli - the list of spikes
    Returns:
        sharedstrs: the two 9 character strings shared between all spikes, at the beginning and end of each spike
        spikedata: a list of dictionaries for each position of a spike, which point a character to spikes with that character in that position
    """
    #spikedata = [['']*(len(spikeli[0][0][1])-18),['']*(len(spikeli[0][0][1])-18)] # initialize - these empty strings will later be dicts
    sharedstrs=[]
    sharedstrs.append((spikeli[0][1][1][:9],spikeli[0][1][1][-9:]))
    sharedstrs.append((spikeli[1][1][1][:9],spikeli[1][1][1][-9:]))
#    print sharedstrs
    spikedata = ''
    #for charno in range(len(spikeli[0][1])-18):
    #    spikedata[charno]=defaultdict(set)
    #    for spike in spikeli:
    #        spikedata[charno][spike[1][charno+9]].add(spike) # build the spikedata dictionary for the specified position
    return sharedstrs,spikedata
    
def conjugate(spikechars):
    newspike = ''
    conjugation = {"A":"T","T":"A","G":"C","C":"G"}
    for char in spikechars[::-1]:
        if char in conjugation:
            newspike += conjugation[char]
        else:
            newspike += char
    return newspike
    
def dist_write(spikedict,path):
    """
    dist_write: writes the spike counts to a file
    Arguments:
        spikedict - the dictionary containing spike counts
        path - the path of the output file
    Returns: nothing
    Effects: creates an output file containing spike counts
    """
    with open(path, 'w') as outfile:
        kvpairs = []
        for key in spikedict:
            if spikedict[key]=="COUNT":
                outfile.write(','.join(list(key)+["COUNT"])+'\n')
            else:
                kvpairs.append((key, spikedict[key]))
        kvpairs = sorted(kvpairs, key=lambda pair: (-pair[1],pair[0][1])) # sort spikes in descending order by count
        for pair in kvpairs:
            outline=list(pair[0])+[str(pair[1])]
            outfile.write(','.join(outline)+'\n') # write the sorted list into the file
    
def main():
    """
    main: sets up command line arguments, then goes through the origin directory to find all fastq files, then calls parse_reads and then dist_write on the file and what is returned from parse_reads, respectively
    Arguments: none, but some arguments do come from the command line (see parser.add_argument help text)
    Returns: nothing
    Effects: Adds the output files to the origin directory
    """
    parser = argparse.ArgumentParser(description='Get inputs for the FASTQ counting script.') # set arguments
    parser.add_argument('spikes', help = 'The file containing spike data.')
    parser.add_argument('source', help = 'The directory to be searched.')
    parser.add_argument('-e','--max_errs',default=0,type=int,help='The maximum number of errors for a spike.')
    parser.add_argument('-r','--removal',action='store_true',help='If specified, removes reads with spikes.')
    args=parser.parse_args()
    
    spikeli = parse_spikes(args.spikes)
    sharedstrs,spikedata=process_spikes(spikeli) # get the spikes and do some preprocessing
    
    for folder, subfolder, filelist in walk(args.source): # look through source folder for fastq files
        for file in filelist:
            name, exten = path.splitext(path.basename(file))
            if exten.lower() == '.fastq': # processes only fastq files
                if folder == args.source:
                    name = folder+name
                else:
                    name = folder + path.sep + name
                spikedict = parse_reads(name, sharedstrs,spikedata,spikeli,args.max_errs,args.removal) # get the spike counts
		#print "spikedict:"
                print 'Processed ' + name
                dist_write(spikedict, name+'.spike.counts.txt') # write the spike counts
    
if __name__=='__main__':
    main()
