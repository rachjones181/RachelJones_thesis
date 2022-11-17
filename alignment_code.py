from msilib import sequence
from operator import length_hint
from Bio import SeqIO
from Bio import pairwise2
import re
from Bio.pairwise2 import format_alignment
from collections import defaultdict  #need to add this to make a dictionary for the mismatches
import pprint


#Force a bad alignment.
#Forces a bad alignment by having a very expensive gap penalty



#starting with all post-degredation samples to get an idea of mutation distribution
#insert sequence to align to here


file locations:
tabv_wt : Projects/RNA_Seq_Proc/data/04_data/fastp_proc/RJ1
tabv+c : Projects/RNA_Seq_Proc/data/03_data/fastp_proc/RJ8/TABV+C_trimmed_alignment.fasta
zika_wt : Projects/RNA_Seq_Proc/data/04_data/fastp_proc/RJ3
zika-c : Projects/RNA_Seq_Proc/data/03_data/fastp_proc/RJ6

tabv_seq  = "GGCAAGGTACGGATTAGCCGTAGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"
tabvc_seq = "GGCAAGGTACGGATTAGCCGTACGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"
                
zika_seq  =  "TGTCAGGCCTGCTAGTCAGCCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"
zikac_seq =  "TGTCAGGCCTGCTAGTCAGCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"

/Users/racheljones/Projects/RNA_Seq_Proc/data/04_data/fastp_proc/RJ3/combined_zwt.fasta

mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RJ1/test_for_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): #this is how you parse thru seqs using SeqIO
        sequence = record.seq #of the many things we can get from SeqIO, just get the actual sequence
        alignments = pairwise2.align.localms(tabv_seq, sequence, 1, -1, -10, -10) #preform alignment, forcing a bad fit that doesn't open gaps
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
                key2 = (one_alignment.seqA[i] + str(i+1))
                value2 = one_alignment.seqB[i]
                mutations_per_position[key2].append(value2)
            else:
                pass


mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RJ2/test_for_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): #this is how you parse thru seqs using SeqIO
        sequence = record.seq #of the many things we can get from SeqIO, just get the actual sequence
        alignments = pairwise2.align.localms(tabvc_seq, sequence, 1, -1, -10, -10) #preform alignment, forcing a bad fit that doesn't open gaps
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
                key2 = (one_alignment.seqA[i] + str(i+1))
                value2 = one_alignment.seqB[i]
                mutations_per_position[key2].append(value2)
            else:
                pass

mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RJ3/test_for_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): #this is how you parse thru seqs using SeqIO
        sequence = record.seq #of the many things we can get from SeqIO, just get the actual sequence
        alignments = pairwise2.align.localms(zika_seq, sequence, 1, -1, -10, -10) #preform alignment, forcing a bad fit that doesn't open gaps
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
                key2 = (one_alignment.seqA[i] + str(i+1))
                value2 = one_alignment.seqB[i]
                mutations_per_position[key2].append(value2)
            else:
                pass

mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RJ4/test_for_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): #this is how you parse thru seqs using SeqIO
        sequence = record.seq #of the many things we can get from SeqIO, just get the actual sequence
        alignments = pairwise2.align.localms(zikac_seq, sequence, 1, -1, -10, -10) #preform alignment, forcing a bad fit that doesn't open gaps
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
                key2 = (one_alignment.seqA[i] + str(i+1))
                value2 = one_alignment.seqB[i]
                mutations_per_position[key2].append(value2)
            else:
                pass

mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RJ4/test_for_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): #this is how you parse thru seqs using SeqIO
        sequence = record.seq #of the many things we can get from SeqIO, just get the actual sequence
        alignments = pairwise2.align.localms(zika_seq, sequence, 1, -1, -10, -10) #perform alignment, forcing a bad fit that doesn't open gaps
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
                key2 = (one_alignment.seqA[i] + str(i+1))
                value2 = one_alignment.seqB[i]
                mutations_per_position[key2].append(value2)
            else:
                pass

#check how our dictionary looks
pprint.pprint(alignment_mismatches)

print(len(alignment_mismatches))
print((alignment_mismatches.values())) #cant count this because it is lists not tuples

#instead iterate through the dictionary and count the number of mismatches per sequence
mutations_per_read = []
for key, value in alignment_mismatches.items():
    mutations_per_read.append(len(value))

for key, value in alignment_mismatches.items():
    if len(value)>20:
        print(key,value)
#do this in terminal so you can copy-paste the numbers into Prism
for i in mutations_per_read:
    print(i)

#####################FIND THE MUTS AT EACH POS###############################
tabv_seq  = "GGCAAGGTACGGATTAGCCGTAGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"
tabvc_seq = "GGCAAGGTACGGATTAGCCGTACGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"

zika_seq  =  "TGTCAGGCCTGCTAGTCAGCCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"
zikac_seq =  "TGTCAGGCCTGCTAGTCAGCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"

#now I need to map what the distribution of mutations actually are at each position:
pprint.pprint(mutations_per_position)
#since it will pull keys completely out of orer, I need to be get them in order so I can copy/paste the data directly
for i in range(0, len(tabv_seq)):
    pos = (tabv_seq[i] + str(i+1))
    #print(pos)
    all_muts = mutations_per_position.get(pos)
    if all_muts is None:
        print("None")
    else: 
        print(all_muts.count("A"))#iterate through all 4 bases manually
        
for i in range(0, len(tabvc_seq)):
    pos = (tabvc_seq[i] + str(i+1))
    #print(pos)
    all_muts = mutations_per_position.get(pos)
    if all_muts is None:
        print("None")
    else: 
        print(all_muts.count("C"))#iterate through all 4 bases manually

for i in range(0, len(zika_seq)):
    pos = (zika_seq[i] + str(i+1))
    #print(pos)
    all_muts = mutations_per_position.get(pos)
    if all_muts is None:
        print("None")
    else: 
        print(all_muts.count("C"))
        
for i in range(0, len(zikac_seq)):
    pos = (zikac_seq[i] + str(i+1))
    #print(pos)
    all_muts = mutations_per_position.get(pos)
    if all_muts is None:
        print("None")
    else: 
        print(all_muts.count("C"))        

print(mutations_per_position.get('C23'))
#this gets put into excel to count totalt muts and frequency
#Now lets say I want to view keys that have a particular mutation I'm interested in and what mutations co-occur with that

for key, value in alignment_mismatches.items():
    if "C23T" in value:
        print(key, value)

#I want to pull out the keys that only had ONE value    
for key, value in alignment_mismatches.items():
    if len(value) < 2 :
        print(value)

#Now I need to make a new dictionary that makes positions as the keys and the numbre of times that mposition gets mutated as the value

for key, value in alignment_mismatches.item():
















########## other code I have modified ###################

alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("RJ3/ZIKA_WT_alignment.fasta") as handle: #read through our trimmed fasta file (all the same length)
    for record in SeqIO.parse(handle, "fasta"): 
        sequence = record.seq
        alignments = pairwise2.align.localms(zika_seq, sequence, 1, -1, -10, -10)
        one_alignment = alignments[0]
        for i in range(0, len(one_alignment.seqA)):
            if one_alignment.seqA[i] != one_alignment.seqB[i]:
                value = ((one_alignment.seqA[i] + str(i+1) + one_alignment.seqB[i]))
                alignment_mismatches[sequence].append(value)
            else:
                pass