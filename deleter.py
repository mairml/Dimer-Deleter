from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import SeqIO
import os
import sys
import regex


adapters=[]
def readadpters(file):
    for primer in SeqIO.parse(file, "fasta"):
        adapter = str(primer.seq)
        adapters.append(adapter)
output = []
def fiveprime(sequences):
    kmer = 2  # number of times to clean the sequence of adapters
    leading = 3  # allowance for bases before adapter appears
    threshold = 3  # allowance for mismatches in the full adapter
    length = 10  # remove adapter if it is at least this length
    for record in SeqIO.parse(sequences, "fastq"):
        sequence = str(record.seq)
        for k in range(kmer,0,-1): #continue doing this for how many possible kmers
            for adapter in adapters:
                for i in range (len(adapter),length-1,-1): #sliding window for degenerate adapters
                    #get the starting region where adapter might appear
                    if k==kmer:
                        string = sequence[0:len(adapter) + leading] #allow leading allowance at start of sequence
                    else:
                        string = sequence[0:len(adapter)]

                    #get end of adapter of various lengths
                    query=adapter[-i:]

                    #organize mismatch allowance
                    if i==length: #only allow one mismatch at smallest size
                        threshold=1
                    threshold=str(threshold)

                    #get the format for regex search
                    matcher = '(?e)('+query+'){e<='+threshold+'}'

                    #regex search
                    match = regex.findall(matcher,string)
                    if len(match)>0: #if the adapter appeared in region
                        indices = [base.start() for base in regex.finditer(matcher, string)]
                        if (indices[0] <= leading) and len(match[0]) == len(query): #get when the adapter first appeared
                            #  #allows for no trailing mismatches
                            qual = record.letter_annotations["phred_quality"]
                            qual = qual[indices[0]+len(query):]
                            sequence = sequence[indices[0]+len(query):] #cut off adapter from sequence
                            record.letter_annotations = {}
                            record.seq = Seq(sequence, generic_dna)
                            record.letter_annotations["phred_quality"] = qual
                            break

        if len(sequence)>=30:
            output.append(record)
    return output


outputF = []
def threeprime(sequences):
    kmer = 1  # number of times to clean the sequence of adapters
    leading = 1  # allowance for bases before adapter appears
    threshold = 3  # allowance for mismatches in the full adapter
    length = 10  # remove adapter if it is at least this length
    for record in SeqIO.parse(sequences, "fastq"):
        sequence = str(record.seq)
        print ("s= "+sequence )
        for k in range(kmer,0,-1): #continue doing this for how many possible kmers
            for adapter in adapters:
                threshold = str(threshold)
                matcher = '(?e)(' + adapter + '){e<=' + threshold + '}'
                #get full match
                indices = [base.start() for base in regex.finditer(matcher, sequence)]
                if len(indices)>0:
                    qual = record.letter_annotations["phred_quality"]
                    qual = qual[0:indices[0]]
                    sequence = sequence[0:indices[0]]  # cut off adapter from sequence
                    record.letter_annotations = {}
                    record.seq = Seq(sequence, generic_dna)
                    record.letter_annotations["phred_quality"] = qual
                    break
                else:
                    #get partial match at end of sequence
                    for i in range (len(adapter)-1,length-1,-1): #sliding window for degenerate adapters
                        #get start of adapter of various lengths
                        query=adapter[0:i]
                        # get the end region where adapter might appear
                        string = sequence[-len(query):]
                        #organize mismatch allowance
                        threshold = 3
                        if i==length: #only allow one mismatch at smallest size
                            threshold=1

                        threshold=str(threshold)

                        #get the format for regex search
                        matcher = '(?e)('+query+'){e<='+threshold+'}'
                        #regex search
                        match = regex.findall(matcher,string)

                        if len(match)>0: #if the adapter appeared in region
                            indices = [base.start() for base in regex.finditer(matcher, sequence)]
                            print ("q= " + query)
                            print ("e= " + string)
                            qual = record.letter_annotations["phred_quality"]
                            qual = qual[0:indices[0]]
                            sequence = sequence[0:indices[0]] #cut off adapter from sequence
                            record.letter_annotations = {}
                            print sequence
                            print indices
                            record.seq = Seq(sequence, generic_dna)
                            record.letter_annotations["phred_quality"] = qual
                            break

        if len(sequence)>=30:
            outputF.append(record)
    return outputF

path = sys.argv[1]
adapterF = sys.argv[2]
out = sys.argv[3]
for f in os.listdir(path):
    file = str(f)
    pfile = path+file
    out = "clean_"+file
    output = [] #five prime
    outputF = [] #three prime
    if ".1" in file:
        print("Trimming sequences of input adapters at the 5' end of "+ file +". Please wait.")
        fiveprime(pfile, adapterF)
        print("Preparing file for output.")
        with open(out, "w") as output_handle:
            SeqIO.write(output, output_handle, "fastq")
    elif ".2" in file:
        print("Trimming 3' end from " + file + ". Please be patient.")
        threeprime(output)
        print("Preparing file for output.")
        with open(out, "w") as output_handle:
            SeqIO.write(outputF, output_handle, "fastq")
print("Program complete!")