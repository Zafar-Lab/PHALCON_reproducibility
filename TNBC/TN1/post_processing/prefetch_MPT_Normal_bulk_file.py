
'''
import os
import pandas as pd
import numpy as np
import re


avinput_dir = '/home/priya/Downloads/Final_Stage/annotation by annovar/Bam post processed files FINAL/Clonal and non clonal distinguished/TN1_bam_processed/Non clonal/TN1_post_bam_nonclonal.avinput'
avinput_df = pd.read_csv(avinput_dir, sep='\t',header=None)

avinput_df.drop([3,4,5,6,7], axis=1, inplace = True)
avinput_df[1] = avinput_df[1] - 1

avinput_df.to_csv('TN1_post_bam_nonclonal.bed',sep='\t',index=None,header=False)

sample_name = 'TN1'
accession_no = 'SRR17032812'
command1 = "prefetch " + accession_no
#os.system(command1)
print("Command finished :", command1)
command2 = "fasterq-dump " + accession_no
#os.system(command2)
print("Command finished :", command2)
command3 = "bowtie2 --no-unal -p 4 -x hg19 -1 "+ accession_no +"_1.fastq -2 " + accession_no + "_2.fastq -S "+ sample_name + ".sam"
os.system(command3)
print("Command finished :", command3)
command4 = "samtools view -@ 4 -Sb -o "+ sample_name +".bam " + sample_name + ".sam"
os.system(command4)
print("Command finished :", command4)
command5 = "samtools sort " + sample_name + ".bam -o " + sample_name + "_sorted.bam"
os.system(command5)
print("Command finished :", command5)
command6 = "samtools mpileup -l " + sample_name + ".bed " + sample_name + "_sorted.bam -o " + sample_name + ".mpileup"
os.system(command6)
print("Command finished :", command6)
command7 = "rm "+ sample_name + ".bam"
os.system(command7)
print("Command finished :", command7)
command8 = "rm "+ sample_name + "_sorted.bam"
os.system(command8)
print("Command finished :", command8)
command9 = "rm "+ accession_no +"_1.fastq"
os.system(command9)
print("Command finished :", command9)
command10 = "rm "+ accession_no +"_2.fastq"
os.system(command10)
print("Command finished :", command10)
command11 = "rm -rf "+ accession_no
os.system(command11)
print("Command finished :", command11)
command12 = "rm "+ sample_name + ".sam"
os.system(command12)
print("Command finished :", command12)

mpileup = pd.read_csv(sample_name + ".mpileup", sep='\t',header=None)
bed_file = pd.read_csv(sample_name + ".bed", sep = '\t', header=None)

mpileup_sites = list(mpileup[1])
print("mpileup sites :",mpileup_sites)



    '''        

import os
import pandas as pd
import numpy as np
import re


avinput_dir = '/home/priya/Downloads/Bowtie/TN5_dbsnp_nz.avinput'
avinput_df = pd.read_csv(avinput_dir, sep='\t',header=None)

avinput_df.drop([3,4,5,6,7], axis=1, inplace = True)
avinput_df[1] = avinput_df[1] - 1

avinput_df.to_csv('TN5_dbsnp_nz.bed',sep='\t',index=None,header=False)

command = "samtools mpileup -l TN5_dbsnp_nz.bed TN5_N_sorted.bam -o TN5_N.mpileup"
os.system(command)

mpileup = pd.read_csv("TN5_N.mpileup", sep='\t',header=None)

mpileup_sites = list(mpileup[1])
print("Tn5 mpileup sites :",mpileup_sites)
