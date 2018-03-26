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
		record_dict= SeqIO.index (reference_genome, "fasta")
		

	except IOError:
		print("Problem opening file")
		exit()
		
	
	READdatafile=generatereads(record_dict,readlength,no_of_reads,error_rate)
	return READdatafile
		
		


def generatereads(record_dict,readlength,no_of_reads,error_rate):
	
	readsfile = open("Readdata1.fasta","w")
	for read in range(no_of_reads):
		readkey=random.choice(list(record_dict))
		
		sequence=[]
		quals='I'* readlength
		sequence = str(record_dict[readkey].seq).upper()
		position_dict={}
		position = random.choice(range(len(sequence)))
		mutated_read=""
		#print(sequence[position+readlength])
		if position + readlength <= (len(sequence)):
			randomread = sequence[position:position+readlength]
			mutated_read=mutate(randomread,error_rate)

		else:
			randomread = sequence[position-readlength:position]
			mutated_read=mutate(randomread,error_rate)

		readsfile.write("@"+readkey)
		readsfile.write("\n")
		readsfile.write(mutated_read)
		readsfile.write("\n")
		readsfile.write("+\n")
		readsfile.write(quals)
		readsfile.write("\n")

	readsfile.close()
	return readsfile
			

def mutate(randomread,error_rate):
	oldread = randomread
	mutatedread=""
	#bases=['A','T','G','C']
	bases='ATGCN'
	for i in range(len(oldread)):
		if random.uniform(0,1) <= error_rate:
			randombase=random.choice(bases)
			if randombase != oldread[i]:
				mutatedread=mutatedread + randombase
			else:
				mutatedread = mutatedread + oldread[i]
		else:
			mutatedread = mutatedread + oldread[i]
	return mutatedread
	

	
		
	
		
## Get the required input from the user like path to reference genome file and readlength

#reference_genome=input("Input the complete file path of reference genome in FASTA format\n")


reference_genome='C:\Pieriandx_assignment\GRCh38_latest_genomic.fna\GRCh38_latest_genomic.fna'


readlength = int(input("Print the readlength for DNA reads\n"))
error_rate = float(input("Enter the desired error rate eg 0.01\n"))
no_of_reads=int(input("Print the no.of reads to generate"))


INPUT_READS=ReadSimulator(reference_genome,readlength,no_of_reads,error_rate)
Mutated_reads = INPUT_READS.name
print(Mutated_reads)
#PerformBWA(reference_genome,Mutated_reads)	