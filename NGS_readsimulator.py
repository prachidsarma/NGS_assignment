#!usr/bin/python
import sys
import Bio
import random
import subprocess
import os

def ReadSimulator(reference_genome,readlength,totalreads,error_rate):
	reference_genome.rstrip("\n")
	try:
		from Bio import SeqIO
		
		record_dict={}
		for seq_record in SeqIO.parse(reference_genome, "fasta"):
			record_dict.update({seq_record.id : seq_record.seq})
		#record_dict = SeqIO.parse(reference_genome,"fasta")
		

	except IOError:
		print("Problem opening file")
		
	
	READdatafile=generatereads(record_dict,readlength,no_of_reads,error_rate)
	return READdatafile
		
		


def generatereads(record_dict,readlength,no_of_reads,error_rate):
	
	readsfile = open("Readdata.fasta","w")
	for read in range(no_of_reads):
		readkey=random.choice(list(record_dict))
		
		sequence=[]
		quals='I'* readlength
		sequence = str(record_dict[readkey]).upper()
		position_dict={}
		position = random.choice(range(len(sequence)))
		
		#print(sequence[position+readlength])
		if position + readlength <= (len(sequence)-1):
			randomread = sequence[position:position+readlength]
			mutated_read=mutate(randomread,error_rate)
			readsfile.write("@"+readkey)
			readsfile.write("\n")
			readsfile.write(mutated_read)
			readsfile.write("\n")
			readsfile.write("+\n")
			readsfile.write(quals)
			readsfile.write("\n")
			position_dict[readkey]=position
	readsfile.close()
	return readsfile
			

def mutate(randomread,error_rate):
	oldread = randomread
	mutatedread=""
	#bases=['A','T','G','C']
	bases='ATGC'
	for i in range(len(oldread)):
		if random.uniform(0,1) < error_rate:
			randombase=random.choice(bases)
			if randombase != oldread[i] and oldread[i] != 'N':
				mutatedread=mutatedread + randombase
				
		else:
			mutatedread = mutatedread + oldread[i]
	return mutatedread
	
def PerformBWA(reference_genome,IN_READS):
	from subprocess import call
	Refgenome_file=reference_genome
	print(reference_genome)
	Input_reads= IN_READS
	#OUTPUT_file = 'Alignedreads.sam'
	#f=open(OUTPUT_file,"w")
	output_sai_file = 'read_1.sai'
	from Bio.Sequencing.Applications import BwaSamseCommandline
	from Bio.Sequencing.Applications import BwaAlignCommandline
	from Bio.Sequencing.Applications import BwaIndexCommandline
	align_cmd=  BwaAlignCommandline(reference=Refgenome_file, read_file=Input_reads)
	cmd = "bwa index reference_genome"

	
		
	
		
## Get the required input from the user like path to reference genome file and readlength

reference_genome=input("Input the complete file path of reference genome in FASTA format\n")


#reference_genome='/home/prachi/Desktop/Prachi/GCF_000001405.25_GRCh37.p13_genomic.fna'


readlength = int(input("Print the readlength for DNA reads\n"))
error_rate = float(input("Enter the desired error rate eg 0.01\n"))
no_of_reads=int(input("Print the no.of reads to generate"))


INPUT_READS=ReadSimulator(reference_genome,readlength,no_of_reads,error_rate)
Mutated_reads = INPUT_READS.name
print(Mutated_reads)
#PerformBWA(reference_genome,Mutated_reads)	