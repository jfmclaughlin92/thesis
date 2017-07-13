#!/usr/bin/env python
# encoding: utf-8

"""
File: SNPmapper.py
Author: Jessica McLaughlin
Last update: 26 December 2016
Info: get coordinates for SNPs in a vcf file relative to location of the UCE probe on the locus
""" 

from Bio import SeqIO
import vcf
import argparse

def get_args():
	"""Get arguments from user"""
	parser = argparse.ArgumentParser(
		description="get coordinates for SNPs in a vcf file relative to location of the UCE probe on the locus"
		)
	parser.add_argument(
		'--input_fasta',
		required=True,
		help='The fasta file to be searched'
	)
	parser.add_argument(
		'--input_vcf',
		required=True,
		help='the input BLASTn xml files'
	)
	parser.add_argument(
		'--probes_list',
		required=True,
		help="the FASTA file with probe sequences"
	)
	parser.add_argument(
		'--output',
		required=True,
		help='the filename for output SNP list'
	)
	return parser.parse_args()
	
def get_probes(probe_file):
	#get probe sequences from file in more usable format
	probes_list=[]
	probe_count=0 # might as well count them too, to make sure everything's working
	
	with open(probe_file,"rU") as probes:
		for probe in SeqIO.parse(probes, "fasta"):
			probes_list.append(probe.seq)
			probe_count+=1
	
	print("Using "+str(probe_count)+" probes.")
	return probes_list	
	
def windows(probe,win_size,step_size):
	# make a sliding window for matching probes to contigs
	probe_len = len(probe)
	
	for i in range(0,probe_len,step_size):
		j = probe_len if i+win_size>probe_len else i+win_size
		yield probe[i:j]
		if j==probe_len: break
		
def find_offset(probe_window,seq,window_loc):
	
	probe_loc = seq.find(probe_window)+window_loc+60 
	#offset is distance from beginning of locus + how many window shifts have occurred + half of probe length (120)
	
	return probe_loc
	
	
def find_probe_locations(probes_list,refseq):

	# match each probe to a contig
	match_count =0 # and count it. just count everything.
	matches=[] # this will be changed to something more sophisticated w/ position info
	with open(refseq, "rU") as ref_fasta:
		for contig in SeqIO.parse(ref_fasta,"fasta"):
			print contig.id
			probe_match=0 # set counter to limit how many records can match
			
			for probe in probes_list:
			 	#print("hi!")
				
				window_loc=0
				for window in windows(probe,50,10):
					if window in contig.seq:
						if probe_match==0:
							offset=find_offset(window,contig.seq,window_loc)
							print offset
							matches.append([contig.id,offset])
							match_count+=1
							probe_match+=1
						else:
							pass	
								
					else:
						window_loc+=10
				
				
						
	print("Found "+str(match_count)+" probe matches in "+str(refseq))
	return matches


def get_SNP_coordinates(locus_list,input_vcf,output):
	try:
		reader = vcf.Reader(open(input_vcf, 'r'))
		print("Reading vcf....")
	except:
			print("Could not open vcf file")
			raise IOError("vcf file "+input_vcf+" could not be opened")
	
	output_file = open(output,'w')
	output_file.write("Locus\tSNP_position\r")
	print("Writing output....")
			
	for record in reader:
		for locus in locus_list:
			if locus[0] in record.CHROM:
			
				new_pos = (record.POS - locus[1])
				output_file.write(str(record.CHROM)+"\t"+str(new_pos)+"\r")
				
	return output_file
	
def main():
	# get args
	args = get_args()

	probe_list= get_probes(args.probes_list) #get probes
	probe_offsets = find_probe_locations(probe_list,args.input_fasta)
	output = get_SNP_coordinates(probe_offsets,args.input_vcf,args.output)

if __name__ == '__main__':
	main()	
	