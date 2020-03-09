#!/usr/bin/python3
import sys
import re
import time
from math import log
from scipy import stats
import numpy as np
import optparse
import configparser
import sys
import os
import utils

parser = optparse.OptionParser(usage = "%prog [options]",version="%prog 1.1.0",description = "Dectect fusion from last-maf")
#parser.add_option("-v","--version",action="store_true",dest="version", default=False,help="print the version")
parser.add_option("-i","--maf",action = "store",type = "string",dest = "maf",default = False,help = "Maf file outputed from last.")
parser.add_option("-o","--output",action = "store",type = "string",dest = "out",default = False,help = "Output file.")
parser.add_option("-c","--config",action = "store",type = "string",dest = "config",default = os.path.dirname(os.path.abspath(__file__))+"/config.ini",help = "Give the full path to your own ini file [default: config.ini].")
(options, args) = parser.parse_args()
if not (options.maf and options.out):
	parser.print_help()
	sys.exit(0)
opt_maf=options.maf
opt_out=options.out
cfg = configparser.ConfigParser()
cfg.read(options.config)
#####get environmental variant
mq_cutoff=int(cfg.get('Filter options', 'mq_cutoff'))
start_cutoff=float(cfg.get('Filter options', 'start_cutoff'))
read_cutoff=int(cfg.get('Detection options', 'read_cutoff'))
merge_bp_cutoff=int(cfg.get('Detection options', 'merge_bp_cutoff'))
merge_cutoff=int(cfg.get('Detection options', 'merge_cutoff'))
precise_cutoff=float(cfg.get('Output filter options', 'precise_cutoff'))
def printINFO(L):
	print("{")
	for i in L.keys():
		if isinstance(i,str) or isinstance(i,int):
			print("%s%s: {" % i)
		else:
			print("%s%s: {" % i.all)
		
		print 
	print("}")
		
def main():
	#####run
	sys.stderr.write(time.strftime("%c") + " Busy with parsing maf file...\n")
	utils.parse_maf.parse_maf()
#	print(">>>read:{")
#	for i in utils.parse_maf.read.keys():
#		seq=""
#		for j in utils.parse_maf.read[i].keys():
#			seq+=" %s: {%s}; " % (j, utils.parse_maf.read[i][j].all)
#		print("%s: {%s}," % (i, seq))
#	print("}\n>>>sv_1:%s" % utils.parse_maf.sv_1)
#	print(">>>sv_2:%s" % utils.parse_maf.sv_2)
	
	sys.stderr.write(time.strftime("%c") + " Busy with print vcf header...\n")
	utils.output_vcf.print_vcf_header()
	
	sys.stderr.write(time.strftime("%c") + " Busy with parsing breakpoints...\n")
	utils.parse_breakpoints.parse_breakpoints()
#	print(">>>key_region:%s" % utils.parse_breakpoints.key_region)
#	print(">>>sv_3:%s" % utils.parse_breakpoints.sv_3)

	sys.stderr.write(time.strftime("%c") + " Busy with parsing svs...\n")
	utils.output_vcf.parse_svs()
	
	
	
if __name__ == "__main__":
	main()
