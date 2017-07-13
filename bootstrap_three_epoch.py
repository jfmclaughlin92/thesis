#!/usr/bin/env python
# encoding: utf-8

"""
File: bootstrap_three_epoch.py
Author: JF McLaughlin
Version: 1.0
Updated: 19 Jan 2017
Info: Given a .vcf file, create multiple .fs input files sampling different individuals and run in DaDi. 
Requirements: convert_vcf_to_dadi_input.pl by Kun Wang, phyluce 1.5, BioPython, PyVCF. And of course, DaDi.
"""

import os
import sys
import shlex
import argparse
import pylab
import dadi
import Demographics1D
import vcf
import collections
import ConfigParser
from Bio import SeqIO
from random import choice
import subprocess 
from numpy import array, random

def get_args():
	"""Get arguments from user"""
	parser = argparse.ArgumentParser(
		description="Given an input vcf file, produce datasets for bootstrapping in dadi"
		)
	parser.add_argument(
		'--input_vcf',
		required=True,
		#action=FullPaths,
		#type=is_file,
		help='The vcf file to be used for creating fs file'
	)
	parser.add_argument(
		'--bootstrap_reps',
		required=True,
		type=int,
		help='How many times to sample data. Greater than 30 recommended.'
	)
	parser.add_argument(
		'--model_conf',
		required=True,
		help='Config file for population info and model parameters'
	)
	parser.add_argument(
		'--output',
		required=True,
		type=str,
		help='Name for output files. Not a full file name, just the base name.'
	)
	parser.add_argument(
		'--indiv_config',
		#type=is_file,
		#action=FullPaths,
		required=True,
		help='List of individual names as appearing in vcf header'
	)	
	parser.add_argument(
		'--reference',
		#type=is_file,
		#action=FullPaths,
		required=True,
		help='fasta file with reference sequence'
	)
	parser.add_argument(
		'--pl_script',
		#type=is_file,
		#action=FullPaths,
		required=True,
		help='where the perl script is'
	)		
	parser.add_argument(
		'--folded',
		#type=is_file,
		#action=FullPaths,
		#required=True,
		default= "YES",
		choices= ["YES","NO"],
		help='folded or unfolded spectra?'
	)
	#parser.add_argument(
		#'--log-path',
		#action=FullPaths,
		#type=is_dir,
		#default=None,
		#help='Directory in which to store log file'
	#)
	return parser.parse_args()
	
def ConfigSectionMap(Config,section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
    
def indiv_by_pop(pop_name,indiv_list):	# get individual names for each population
	taxon=[] # set up empty list
	for ind in indiv_list:
		if pop_name in ind:
			taxon.append(ind) 
	return taxon
    
def data_to_lines(data):
	return data.split("\n")
	
def bootstrap_sample(indiv_list, pop_size, i_id):

	bootstrap_indiv = random.choice(indiv_list, pop_size)
	print("Bootstrap "+str(i_id)+"using: "+str(bootstrap_indiv))
	
	return bootstrap_indiv

def parse_vcf(vcf_input,output_tag, individuals, i_id):
	try:
		reader = vcf.Reader(open(vcf_input, 'r'))
		print("Reading vcf....")
	except:
			print("Could not open vcf file")
			raise IOError("vcf file "+vcf_input+" could not be opened")
	
	clean_records =[]		
	for record in reader: # for each record
		selected_calls = []
		
		for call in record.samples: 
			for indiv in individuals: 	
			# go to our selected individuals
				if indiv in call.sample: # if it's a match...
						selected_calls.append(call)
				else:
						pass 

		record.samples = selected_calls
		clean_records.append(record)
		
	use_samples = []
	for sample in reader.samples:
		for indiv in individuals:
			if indiv in sample:
				use_samples.append(sample)
	reader.samples = use_samples

	try:
		output_vcf = output_tag+'.vcf'
		writer = vcf.Writer(open(output_vcf, 'w'), reader) 
		print('Output vcf created')
	except:
			print('Could not create output vcf file')
			raise IOError('Could not create output vcf file')
				

	for record in clean_records:
		writer.write_record(record)		
	
	return output_vcf
	
def create_fs_list(individuals,out_tag,taxon,pop):
	# create the list file necessary to run the vcf to dadi input script
	list_file = open(out_tag+"list.txt", 'w')
	try:
		for indiv in individuals:
			if indiv in taxon:
				list_file.write(indiv+"\t"+pop+"\n")
			else:
				pass
	except:
		raise IOError("Could not write list file")
		
	list_file.close()
	
	return list_file

			
def create_fs(perl_script,vcf_input,refseq,output_tag,listfile,pop,count,folded):

	args_str = "perl "+perl_script+" "+refseq+" "+vcf_input+" "+listfile.name
	print args_str
	args = shlex.split(args_str)
	print args
#	perl_input = ['convert_vcf_to_dadi_input.pl', refseq, vcf_input, listfile]
#	cmd = [
#		"perl",
#		perl_script,
#		refseq,
#		vcf_input,
#		listfile
#	]
#	print cmd 

	perl_out = subprocess.Popen(args, stdout=subprocess.PIPE)
	perl_out.communicate()
	print("vcf converted to dadi input successfully")
	
	#		print("Could not run .vcf to .fs conversion script") # kill if script not used
	#		raise IOError("Could not create fs-type file")
	

	dd = dadi.Misc.make_data_dict(listfile.name+".data") # use dadi to make a data dictionary for this data
	if folded == "YES":
		folded = dadi.Spectrum.from_data_dict(dd,[pop],[count],polarized=False) # generate a folded spectrum
		fs_file = folded.to_file(output_tag+'.fs') # use this to make a folded file
		
	elif folded == "NO":
		fs = dadi.Spectrum.from_data_dict(dd,[pop],[count],polarized=True) # for unfolded Aegolius owls
	
		fs_file = fs.to_file(output_tag+'.fs') # use this to make an unfolded file	
	print("fs file complete!")
	return fs_file



def run_dadi(fs_file,output_tag,i_id,dadi_iter,nuB,nuF,TB,TF):
	#TO DO: add functions to allow other model types
 	# Load the data
	data = dadi.Spectrum.from_file(output_tag+".fs")
	ns = data.sample_sizes
	data.mask_corners()

	# create a file to print output to and set number of iterations
	output_file = open(output_tag+".txt", "w")
	iterations = dadi_iter

	# print number of samples to verify correct load
	print(ns)

	# Grid point settings will be used for extrapolation.
	# Grid points need to be formated [n,n+10,n+20]. 
	# Needs to be bigger than the number of samples you have (n>ns) and this will be a strong determination as to how long your program will run.

	pts_l = [110,120,130]

	# Call particular model to run, the model choosen here is split.w.migration
	func = dadi.Demographics1D.three_epoch

	# Load the data
	#	
	#   
	#	params = (nuB,nuF,TB,TF)   
	# 	ns = (n1,)
	#
	#	nuB: Ratio of bottleneck population size to ancient pop size
	#	nuF: Ratio of contemporary to ancient pop size
	#	TB: Length of bottleneck (in units of 2*Na generations) 
	#	TF: Time since bottleneck recovery (in units of 2*Na generations) 

	# Now let's optimize parameters for this model.

	# The upper_bound and lower_bound lists are for use in optimization.
	# Occasionally the optimizer will try wacky parameter values. We in particular
	# want to exclude values with very long times, very small population sizes, or
	# very high migration rates, as they will take a long time to evaluate.
	# Parameters are: (nuB,nuF,TB,TF)
	
	#Set the upper and lower bounds to make sure that the bounderies are 
	#there. Suggested time parameters: lower 0, upper 5, migration 
	#parameters: lower 0, upper 10,size parameters: lower 1e-2, upper 100 
	upper_bound = [nuB,nuF,TB,TF]
	lower_bound = [1e-1, 1e-1, 0, 0]

	# ready to run multiple iterations
	output_file.write("Iteration,MLCL,nuB,nuF,TB,TF,theta\n") 
	

	for iter in range(iterations):
		# create an id for this iteration for labelling purposes
		iter_id = str(iter)
	
		# This is our initial guess for the parameters, which is somewhat arbitrary.
		p0 = [1,1,1,1]
		# Make the extrapolating version of our demographic model function.
		func_ex = dadi.Numerics.make_extrap_log_func(func)

		# Perturb our parameters before optimization. This does so by taking each
		# parameter up to a factor of two up or down.
		p0 = dadi.Misc.perturb_params(p0, fold=2, upper_bound=upper_bound,
 	    	  	                      lower_bound=lower_bound)

		print('Beginning optimization ************************************************')
		popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
    		                               lower_bound=lower_bound,
        		                           upper_bound=upper_bound,
            		                       verbose=4)
		# The verbose argument controls how often progress of the optimizer should be
		# printed. It's useful to keep track of optimization process.
		print('Finshed optimization **************************************************')

		print('Best-fit parameters: {0}'.format(popt))
		

		# Calculate the best-fit model AFS.
		model = func_ex(popt, ns, pts_l)

		# Likelihood of the data given the model AFS.
		ll_model = dadi.Inference.ll_multinom(model, data)
		print('Maximum log composite likelihood: {0}'.format(ll_model))

		# The optimal value of theta given the model.
		theta = dadi.Inference.optimal_sfs_scaling(model, data)
		print('Optimal value of theta: {0}'.format(theta))
		print(str(pts_l)+", "+str(upper_bound)+", "+str(lower_bound))
		output_file.write(iter_id+",{0}".format(ll_model)+'{0}'.format(popt)+',{0}'.format(theta)+"\n")
		
	
		# Plot a comparison of the resulting fs with the data.
		
		pylab.figure(1)
		dadi.Plotting.plot_1d_comp_multinom(model, data)

		# This saves the resulting figure to file
		pylab.savefig(output_tag+iter_id+"_"+str(i_id)+'.png')
	output_file.close()
	return output_file
		
def main():
	# get args
	args = get_args()
    # setup logging
	#log, my_name = setup_logging(args) 
	print("Logging set-up")
	
	print("Parsing conf file......")
	Config = ConfigParser.ConfigParser()
	Config.read(args.model_conf)
	print("Getting config info for data........")

	pop = ConfigSectionMap(Config,"Data")["population"]

	N = int(ConfigSectionMap(Config, "Data")["pop_size"])


	n_alleles_1 = N * 2

	pop_alleles = [n_alleles_1]

	print("Getting config info for model.........")
	
	# get the info from the conf file. TO DO: create as separate functions for different model types
	NuB = float(ConfigSectionMap(Config, "Model")["nub"])
	NuF = float(ConfigSectionMap(Config, "Model")["nuf"])
	TB = float(ConfigSectionMap(Config, "Model")["tb"])
	TF = float(ConfigSectionMap(Config, "Model")["tf"])
	dadi_reps = int(ConfigSectionMap(Config, "Model")["dadi_runs"])
	upper_bound = [NuB,NuF,TB,TF]
	
	print("Using upper-bound values "+str(upper_bound))
	print("Dadi reps requested: "+str(dadi_reps))

	
	print("Getting list of unique individuals")
	individual_file = open(args.indiv_config,'r')
	individuals_data = individual_file.read()
	individuals = data_to_lines(individuals_data)
	taxon = indiv_by_pop(pop, individuals)

	folded= args.folded
	
	bootstraps = range(args.bootstrap_reps)
	out_tag = args.output
	print out_tag
	vcf_in = args.input_vcf
	refseq = args.reference
	print refseq
	
	for b in bootstraps:
		print("Bootstrap run: "+str(b))

		if not os.path.isdir(out_tag+str(b)): #make an output directory for this round
			os.makedirs(out_tag+str(b))
		else:
			pass
		out_file = out_tag+str(b)+"/"+out_tag+str(b)
		print out_file
		bs_indiv = bootstrap_sample(taxon, N, b).tolist()
		bs_vcf = parse_vcf(vcf_in, out_file, bs_indiv,b)
		bs_fs_list = create_fs_list(bs_indiv,out_file,taxon,pop)
		bs_fs = create_fs(args.pl_script,bs_vcf,refseq,out_file,bs_fs_list,pop,n_alleles_1,folded)
		dadi_runs = run_dadi(bs_fs,out_file,b,dadi_reps,NuB,NuF,TB,TF)
	
	print("Done!")
	
if __name__ == '__main__':
	main()	
