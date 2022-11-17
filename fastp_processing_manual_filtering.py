activate environment : 
source /Users/racheljones/.virtualenvs/env2/bin/activate
cd /Users/racheljones/Projects/RNA_Seq_Proc/data/


#Initial data process step, connects paired end reads and removes low quality reads (not the correct file names)
>1: RAJ_01
fastp -i raw/RAJ_01/RAJ_01_R1.fq.gz -o fastp_proc/RAJ_01/out.RAJ_01_R1.fq.gz -I raw/RAJ_01/RAJ_01_R2.fq.gz -O fastp_proc/RAJ_01/out.RAJ_01_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_01/RAJ_01_merged.fq.gz #--reads_to_process 5000
>2: RAJ_02
fastp -i raw/RAJ_02/RAJ_02_R1.fq.gz -o fastp_proc/RAJ_02/out.RAJ_02_R1.fq.gz -I raw/RAJ_02/RAJ_02_R2.fq.gz -O fastp_proc/RAJ_02/out.RAJ_02_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_02/RAJ_02_merged.fq.gz #--reads_to_process 5000
>3: RAJ_03
fastp -i raw/RAJ_03/RAJ_03_R1.fq.gz -o fastp_proc/RAJ_03/out.RAJ_03_R1.fq.gz -I raw/RAJ_03/RAJ_03_R2.fq.gz -O fastp_proc/RAJ_03/out.RAJ_03_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_03/RAJ_03_merged.fq.gz #--reads_to_process 5000
>4: RAJ_04
fastp -i raw/RAJ_04/RAJ_04_R1.fq.gz -o fastp_proc/RAJ_04/out.RAJ_04_R1.fq.gz -I raw/RAJ_04/RAJ_04_R2.fq.gz -O fastp_proc/RAJ_04/out.RAJ_04_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_04/RAJ_04_merged.fq.gz #--reads_to_process 5000
>5: RAJ_05
fastp -i raw/RAJ_05/RAJ_05_R1.fq.gz -o fastp_proc/RAJ_05/out.RAJ_05_R1.fq.gz -I raw/RAJ_05/RAJ_05_R2.fq.gz -O fastp_proc/RAJ_05/out.RAJ_05_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_05/RAJ_05_merged.fq.gz #--reads_to_process 5000
>6: RAJ_06
fastp -i raw/RAJ_06/RAJ_06_R1.fq.gz -o fastp_proc/RAJ_06/out.RAJ_06_R1.fq.gz -I raw/RAJ_06/RAJ_06_R2.fq.gz -O fastp_proc/RAJ_06/out.RAJ_06_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_06/RAJ_06_merged.fq.gz #--reads_to_process 5000
>7: RAJ_07
fastp -i raw/RAJ_07/RAJ_07_R1.fq.gz -o fastp_proc/RAJ_07/out.RAJ_07_R1.fq.gz -I raw/RAJ_07/RAJ_07_R2.fq.gz -O fastp_proc/RAJ_07/out.RAJ_07_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_07/RAJ_07_merged.fq.gz #--reads_to_process 5000
>8: RAJ_08
fastp -i raw/RAJ_08/RAJ_08_R1.fq.gz -o fastp_proc/RAJ_08/out.RAJ_08_R1.fq.gz -I raw/RAJ_08/RAJ_08_R2.fq.gz -O fastp_proc/RAJ_08/out.RAJ_08_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_08/RAJ_08_merged.fq.gz #--reads_to_process 5000
>9: RAJ_09
fastp -i raw/RAJ_09/RAJ_09_R1.fq.gz -o fastp_proc/RAJ_09/out.RAJ_09_R1.fq.gz -I raw/RAJ_09/RAJ_09_R2.fq.gz -O fastp_proc/RAJ_09/out.RAJ_09_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_09/RAJ_09_merged.fq.gz #--reads_to_process 5000
>10: RAJ_10
fastp -i raw/RAJ_10/RAJ_10_R1.fq.gz -o fastp_proc/RAJ_10/out.RAJ_10_R1.fq.gz -I raw/RAJ_10/RAJ_10_R2.fq.gz -O fastp_proc/RAJ_10/out.RAJ_10_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_10/RAJ_10_merged.fq.gz #--reads_to_process 5000
>11: RAJ_11
fastp -i raw/RAJ_11/RAJ_11_R1.fq.gz -o fastp_proc/RAJ_11/out.RAJ_11_R1.fq.gz -I raw/RAJ_11/RAJ_11_R2.fq.gz -O fastp_proc/RAJ_11/out.RAJ_11_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_11/RAJ_11_merged.fq.gz #--reads_to_process 5000
>12: RAJ_12
fastp -i raw/RAJ_12/RAJ_12_R1.fq.gz -o fastp_proc/RAJ_12/out.RAJ_12_R1.fq.gz -I raw/RAJ_12/RAJ_12_R2.fq.gz -O fastp_proc/RAJ_12/out.RAJ_12_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_12/RAJ_12_merged.fq.gz #--reads_to_process 5000
>13: RAJ_13
fastp -i raw/RAJ_13/RAJ_13_R1.fq.gz -o fastp_proc/RAJ_13/out.RAJ_13_R1.fq.gz -I raw/RAJ_13/RAJ_13_R2.fq.gz -O fastp_proc/RAJ_13/out.RAJ_13_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_13/RAJ_13_merged.fq.gz #--reads_to_process 5000
>14: RAJ_14
fastp -i raw/RAJ_14/RAJ_14_R1.fq.gz -o fastp_proc/RAJ_14/out.RAJ_14_R1.fq.gz -I raw/RAJ_14/RAJ_14_R2.fq.gz -O fastp_proc/RAJ_14/out.RAJ_14_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_14/RAJ_14_merged.fq.gz #--reads_to_process 5000
>15: RAJ_15
fastp -i raw/RAJ_15/RAJ_15_R1.fq.gz -o fastp_proc/RAJ_15/out.RAJ_15_R1.fq.gz -I raw/RAJ_15/RAJ_15_R2.fq.gz -O fastp_proc/RAJ_15/out.RAJ_15_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_15/RAJ_15_merged.fq.gz #--reads_to_process 5000
>16: RAJ_16
fastp -i raw/RAJ_16/RAJ_16_R1.fq.gz -o fastp_proc/RAJ_16/out.RAJ_16_R1.fq.gz -I raw/RAJ_16/RAJ_16_R2.fq.gz -O fastp_proc/RAJ_16/out.RAJ_16_R2.fq.gz  --stdout --merge --merged_out fastp_proc/RAJ_16/RAJ_16_merged.fq.gz #--reads_to_process 5000


#Next, I'll filter out anything smaller than 100 nucleotides to lower the file size and remove degraded RNA fragments

>1: RAJ_01
fastp -i fastp_proc/RAJ_01/RAJ_01_merged.fq.gz -o fastp_proc/RAJ_01/RAJ_01_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>2: RAJ_02
fastp -i fastp_proc/RAJ_02/RAJ_02_merged.fq.gz -o fastp_proc/RAJ_02/RAJ_02_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>3: RAJ_03
fastp -i fastp_proc/RAJ_03/RAJ_03_merged.fq.gz -o fastp_proc/RAJ_03/RAJ_03_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>4: RAJ_04
fastp -i fastp_proc/RAJ_04/RAJ_04_merged.fq.gz -o fastp_proc/RAJ_04/RAJ_04_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>5: RAJ_05
fastp -i fastp_proc/RAJ_05/RAJ_05_merged.fq.gz -o fastp_proc/RAJ_05/RAJ_05_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>6: RAJ_06
fastp -i fastp_proc/RAJ_06/RAJ_06_merged.fq.gz -o fastp_proc/RAJ_06/RAJ_06_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>7: RAJ_07
fastp -i fastp_proc/RAJ_07/RAJ_07_merged.fq.gz -o fastp_proc/RAJ_07/RAJ_07_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>8: RAJ_08
fastp -i fastp_proc/RAJ_08/RAJ_08_merged.fq.gz -o fastp_proc/RAJ_08/RAJ_08_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>9: RAJ_09
fastp -i fastp_proc/RAJ_09/RAJ_09_merged.fq.gz -o fastp_proc/RAJ_09/RAJ_09_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>10: RAJ_10
fastp -i fastp_proc/RAJ_10/RAJ_10_merged.fq.gz -o fastp_proc/RAJ_10/RAJ_10_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>11: RAJ_11
fastp -i fastp_proc/RAJ_11/RAJ_11_merged.fq.gz -o fastp_proc/RAJ_11/RAJ_11_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>12: RAJ_12
fastp -i fastp_proc/RAJ_12/RAJ_12_merged.fq.gz -o fastp_proc/RAJ_12/RAJ_12_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>13: RAJ_13
fastp -i fastp_proc/RAJ_13/RAJ_13_merged.fq.gz -o fastp_proc/RAJ_13/RAJ_13_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>14: RAJ_14
fastp -i fastp_proc/RAJ_14/RAJ_14_merged.fq.gz -o fastp_proc/RAJ_14/RAJ_14_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>15: RAJ_15
fastp -i fastp_proc/RAJ_15/RAJ_15_merged.fq.gz -o fastp_proc/RAJ_15/RAJ_15_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000
>16: RAJ_16
fastp -i fastp_proc/RAJ_16/RAJ_16_merged.fq.gz -o fastp_proc/RAJ_16/RAJ_16_trimmed_merged.fastq --dedup --length_required 100 #--reads_to_process 5000


#TABV_WT and TABV+C
import collections
from itertools import islice
import sys 

newfile = open("fastp_proc/RAJ_06/RAJ_06_filtered_trimmed_merged.fastq", "w")
with open("fastp_proc/RAJ_06/RAJ_06_trimmed_merged.fastq") as file:
    before = collections.deque(maxlen=1)
    for line in file:
        if "TTAGCTTTG" in line: #part of zika leeader
            continue
        if "AATTTTAA" in line: 
            continue
        if "TAAGAGAA" in line: #If the leader sequence is present, exit the loop
            continue
        if "GTGCTGTA" in line: #zika leader, which accounts for some contamination
            continue
        if "GCACCAA" in line:
            continue
        if "TTGTTAGGTAG" in line:
            continue
        else:
            if "GAGGAA" in line: #if the adapter sequence is present, proceed
                if line.startswith("@"):
                    continue
                #trimmed_seq=(line.split("ATCT")[1]) #trim off the adapter
                n = len(line) #whats the length of the remaining sequence
                if n < 50 or n > 115:  #I was getting issuess with short leftover reads not being written properly so filter the shortest reads
                    continue
                newfile.writelines(before)
                newfile.write(line + "\n")
                newfile.writelines(next(file))
                newfile.writelines(next(file)[-n::])
            else:
                pass
            before.append(line)
    newfile.close()

newfile = open("fastp_proc/RAJ_09/RAJ_09_filtered_trimmed_merged.fastq", "w")
with open("fastp_proc/RAJ_09/RAJ_09_trimmed_merged.fastq") as file:
    before = collections.deque(maxlen=1)
    for line in file:
        if "TTGTTAGGTAG" in line:
            continue
        if "TTTTAAGAG" in line:
            continue
        if "CCAATCTTAATGT" in line: #ZiKV lead
            continue 
        else:
            if "GCTCTTCCGATCT" in line: #if the adapter sequence is present, proceed
                if line.startswith("@"):
                    continue
                trimmed_seq=(line.split("GCTCTTCCGATCT")[1]) #trim off the adapter
                n = len(trimmed_seq) #whats the length of the remaining sequence
                #n = len(line)
                if n < 100 or n > 102 :  #I was getting issuess with short leftover reads not being written properly so filter the shortest reads
                    continue
                newfile.writelines(before)
                newfile.write(trimmed_seq)
                newfile.writelines(next(file))
                newfile.writelines(next(file)[-n::])
            else:
                pass
            before.append(line)
    newfile.close()        


#####TESTING FOR SOMETHING
newfile = open("fastp_proc/RAJ_09/RAJ_09_test.fastq", "w")
with open("fastp_proc/RAJ_09/RAJ_09_trimmed_merged.fastq") as file:
    before = collections.deque(maxlen=1)
    for line in file:
        if "TTGTTAGGTAG" in line:
            continue
        if "TTTTAAGAG" in line:
            continue
        if "CCAATCTTAATGT" in line: #ZiKV lead
            continue 
        else:
            if "TCCGATCT" in line: #if the adapter sequence is present, proceed
                if line.startswith("@"):
                    continue
                #trimmed_seq=(line.split("ATCT")[1]) #trim off the adapter
                n = len(line) #whats the length of the remaining sequence
                if n < 25:  #I was getting issuess with short leftover reads not being written properly so filter the shortest reads
                    continue
                newfile.writelines(before)
                newfile.write(line + "\n")
                newfile.writelines(next(file))
                newfile.writelines(next(file)[-n::])
            else:
                pass
            before.append(line)
    newfile.close()

#that fix made it so that there were now some empy lines tho, so I rewrote the file all over again to remove those.          
newfile = open("fastp_proc/RAJ_06/RAJ_06_filtered_trimmed_merged2.fastq", "w")
with open("fastp_proc/RAJ_06/RAJ_06_filtered_trimmed_merged.fastq") as file:
    for line in file:
        if not line.strip():
            continue
        else:
            newfile.write(line[-105:])
    newfile.close()            
#####do not mess with above this line

>2:
fastaptamer_count -i fastp_proc/RAJ_02/RAJ_02_filtered_trimmed_merged2_fixed.fastq -o fastp_proc/RAJ_02/RAJ_02_fasta.fasta

#ZIKA WT and ZIKA-C

newfile = open("fastp_proc/RAJ_03/RAJ_03_filtered_trimmed_merged.fastq", "w")
with open("fastp_proc/RAJ_03/RAJ_03_trimmed_merged.fastq") as file:
    before = collections.deque(maxlen=1)
    for line in file:
        if "ACCAATCTTAA" in line:
            continue
        else:
            if "TCCGATCT" in line: #if the adapter sequence is present, proceed
                if line.startswith("@"):
                    continue
                #trimmed_seq=(line.split("ATCT")[1]) #trim off the adapter
                n = len(line) #whats the length of the remaining sequence
                if n < 25:  #I was getting issuess with short leftover reads not being written properly so filter the shortest reads
                    continue
                newfile.writelines(before)
                newfile.write(line + "\n")
                newfile.writelines(next(file))
                newfile.writelines(next(file)[-n::])
            else:
                pass
            before.append(line)
    newfile.close()

newfile = open("fastp_proc/RAJ_08/RAJ_08_filtered_trimmed_merged2.fastq", "w")
with open("fastp_proc/RAJ_08/RAJ_08_filtered_trimmed_merged.fastq") as file:
    for line in file:
        if not line.strip():
            continue
        else:
            newfile.write(line)
    newfile.close()



####after fasta make sure all are same length and align with alignment_code######

from msilib import sequence
from operator import length_hint
from Bio import SeqIO
from Bio import pairwise2
import re
from Bio.pairwise2 import format_alignment
from collections import defaultdict  #need to add this to make a dictionary for the mismatches
import pprint


tabv_seq  = "GGCAAGGTACGGATTAGCCGTAGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"
tabvc_seq = "GGCAAGGTACGGATTAGCCGTACGGGGCTTGAGAACCCCCCCTCCCCACTCATTTTATTTCCTC"
                
zika_seq  =  "TGTCAGGCCTGCTAGTCAGCCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"
zikac_seq =  "TGTCAGGCCTGCTAGTCAGCACAGCTTGGGGAAAGCTGTGCAGCCTGTGACCCCCCC"


mutations_per_position = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
alignment_mismatches = defaultdict(list) #begin by creating a defaultsict that will store sequence/mismatch pairs
with open("fastp_proc/RAJ_02/RAJ_02_fasta.fasta") as handle: #read through our trimmed fasta file (all the same length)
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
#check how our dictionary looks
pprint.pprint(alignment_mismatches)

print(len(alignment_mismatches))
print((alignment_mismatches.values())) #cant count this because it is lists not tuples

#instead iterate through the dictionary and count the number of mismatches per sequence
mutations_per_read = []
for key, value in alignment_mismatches.items():
    mutations_per_read.append(len(value))

for key, value in alignment_mismatches.items():
    if len(value)>30:
        print(key,value)
#do this in terminal so you can copy-paste the numbers into Prism
for i in mutations_per_read:
    print(i)
