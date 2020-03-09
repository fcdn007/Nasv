#!/usr/bin/python3
import io
import sys
import vcf as py_vcf
import os
import time
import re
from classes import *
from utils import function as f
from utils.parse_breakpoints import sv_3

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NaSV
def print_vcf_header():
	"""Creates vcf header by setting vcf format and calling functions to create each header section."""
	outf=open(NaSV.opt_out,'w')
	outf.write('##fileformat=VCFv4.2\n##fileDate=%s##source=NaSV\n' % str(time.strftime("%c")))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("DEL", "Deletion"))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("DUP", "Duplication"))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("INS", "Insertion"))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("INV", "Invertion"))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("TRA", "Translocation"))
	outf.write('##ALT=<ID=%s,Description="%s">\n' % ("BND", "Breakend"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ('IMPRECISE', 0, "Flag", "Imprecise structural variant"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("SVTYPE", 1, "String", "Type of structural variant"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("EVENT", 1, "String", "ID of event associated to breakend"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("SVMETHOD", 1, "String", "Type of approach used to detect SV"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("CHR2", 1, "String", "Chromosome for END coordinate in case of a translocation"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("END", 1, "Integer", "End position of structural variant"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("CIPOS", 2, "Integer", "Confidence interval (95%) around POS"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("CIEND", 2, "Integer", "Confidence interval (95%) around END"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("SVLEN", '.', "Integer", "Distance between the two genomic positions"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("INSLEN", '.', "Integer", "Insertion length"))
	outf.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ("STRANDS", '.', "String","Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)"))
	outf.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ('GT', 1, 'String', 'Genotype'))
	outf.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ('SU', 2, 'Integer', 'Number of pieces of evidence supporting the variant'))
	outf.write('##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % ('SR', 2, 'Integer', 'Number of split reads supporting the variant'))
	outf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n')
	outf.close()

def parse_svs():
	outf=open(NaSV.opt_out,'a')
	event=0
	for i in sorted(sv_3.keys(),key=lambda x: len(x.split('|')), reverse=True):
		event+=1
		sv_out={}
		next_=0
		for j in range(1,len(sv_3[i]['region'])):
			if next_:
				continue
			key1=sv_3[i]['region'][j-1].keys()
			key2=sv_3[i]['region'][j].keys()
			pos1=[]
			gap=[]
			pos2=[]
			imprecise,SU,gap_j=(1,0,0)
			for k in key1:
				if not isinstance(k,int):
					l_k=re.split(':|_',k)
					strand1=l_k[2][-1]
					if not re.match('NONE',l_k[2]):
						pos1+=[ChrPos(l_k[0], l_k[2][:-1])] * sv_3[i]['region'][j-1][k] 
				SU+=sv_3[i]['region'][j-1][k]
			for k in key2:
				if isinstance(k, int):
					gap_j+=1
			if gap_j==len(key2):
				next_=1
				for k in key2:
					gap+=[int(k)] * sv_3[i]['region'][j][k]
					SU+=sv_3[i]['region'][j][k]
				key2=sv_3[i]['region'][j+1].keys()
				for k in key2:
					if not isinstance(k,int):
						l_k=re.split(':|_',k)
						strand2=l_k[2][-1]
						if not re.match('NONE',l_k[1]):
							pos2+=[ChrPos(l_k[0], l_k[1])] * sv_3[i]['region'][j+1][k] 
					SU+=sv_3[i]['region'][j+1][k]
				SU/=3
			else:
				for k in key2:
					SU+=sv_3[i]['region'][j][k]
					if not isinstance(k, int):
						l_k=re.split(':|_',k)
						strand2=l_k[2][-1]
						if not re.match('NONE',l_k[1]):
							pos2+=[ChrPos(l_k[0], l_k[1])] * sv_3[i]['region'][j][k] 
				SU/=2
			if strand2=="+":
				strand2="-"
			else:
				strand2="+"
			pos1_real,p1=f.representative(pos1)
			gap_real=0
			if len(gap):
				gap_real,_=f.representative(gap)
			pos2_real,p2=f.representative(pos2)
			CIpos1=[]
			CIpos2=[]
			sv_out[j+1]={}
			sv_out[j+1]['STRAND']="".join([strand1,strand2])
			sv_out[j+1]['pos1']=pos1_real
			sv_out[j+1]['pos2']=pos2_real
			sv_out[j+1]['SR']=sv_3[i]['SWR']
			sv_out[j+1]['SU']=SU
			if p1==p2==1:
				imprecise=0
			else:
				CIpos1=f.CI(pos1)
				CIpos2=f.CI(pos2)               
			sv_out[j+1]['imprecise']=imprecise
#####decide SVTYPE
			if gap_real and gap_real<=NaSV.merge_cutoff*2:#ingore small gap
				pos2_real.pos=pos2_real.pos-gap_real if strand2=='-' else pos2_real.pos+gap_real
				gap_real=0
			elif gap_real and pos1_real.chr==pos2_real.chr and pos2_real.pos-pos1_real.pos<gap_real and strand1=="+" and strand2=='-':#decide ins
				gap_real-=pos2_real.pos-pos1_real.pos
				sv_out[j+1]['SVTYPE']='INS' 
			elif gap_real and pos1_real.chr==pos2_real.chr and pos1_real.pos-pos2_real.pos<gap_real and strand1=="-" and strand2=='+':#decide ins
				gap_real-=pos2_real.pos-pos1_real.pos
				sv_out[j+1]['SVTYPE']='INS'
        
			if gap_real==0 and pos1_real.chr==pos2_real.chr and pos2_real.pos-pos1_real.pos>0 and strand1=="+" and strand2=='-':#decide del
				sv_out[j+1]['SVTYPE']='DEL'
			elif gap_real==0 and pos1_real.chr==pos2_real.chr and pos1_real.pos-pos2_real.pos>0 and strand1=="-" and strand2=='+':#decide del
				sv_out[j+1]['SVTYPE']='DEL'

			if strand1=="+" and strand2=="-" and pos1_real.chr==pos2_real.chr and pos1_real.pos>pos2_real.pos:
				sv_out[j+1]['SVTYPE']='DUP'
			elif strand1=="-" and strand2=="+" and pos1_real.chr==pos2_real.chr and pos1_real.pos<pos2_real.pos:
				sv_out[j+1]['SVTYPE']='DUP'
        
			if 'SVTYPE' not in sv_out[j+1]:
				sv_out[j+1]['SVTYPE']='BND'

			if imprecise:
				sv_out[j+1]['CIpos1']=CIpos1
				sv_out[j+1]['CIpos2']=CIpos2
			if gap_real!=0:
				sv_out[j+1]['gap']=gap_real
		next_=0
		print_idx=1
		for e in sv_out.keys():
			if next_:
				next_=0
				continue
#decide INV
			pos11=sv_out[e]['pos1']
			pos12=sv_out[e]['pos2']
			output=[pos11.chr,pos11.pos,'_'.join([str(event),str(print_idx)]),'N']
			if e+1 in sv_out:
				pos21=sv_out[e+1]['pos1']
				pos22=sv_out[e+1]['pos2']
				if pos11.chr==pos12.chr==pos21.chr==pos22.chr and 'gap' not in sv_out[e] and 'gap' not in sv_out[e+1] \
				   and sv_out[e]['STRAND'][0]!=sv_out[e+1]['STRAND'][1] and sv_out[e]['STRAND'][0]==sv_out[e]['STRAND'][1] \
				   and abs(pos11.pos-pos21.pos)<=NaSV.merge_cutoff and abs(pos22.pos-pos12.pos)<=NaSV.merge_cutoff:
					next_=1
					output+=['<INV>','.','.']
					if sv_out[e]['imprecise']==0:
						output+=['SVTYPE=INV;EVENT=%s;CHR2=%s;END=%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, abs(pos12.pos-pos11.pos))]
					else:
						output+=['IMPRECISE;SVTYPE=INV;EVENT=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s;SVLEN=%s' % (event, pos12.chr, \
						pos12.pos, sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],sv_out[e]['CIpos2'][1], abs(pos12.pos-pos11.pos))]
					output+=['GT:SU:SR','./.:%d:%d' % (sv_out[e]['SU'], sv_out[e]['SR'])]
					outf.write('\t'.join([str(x) for x in output])+'\n')
					print_idx+=1
					continue
			if(sv_out[e]['SVTYPE']=='INS'):
				output+=['<INS>','.','.']
				if sv_out[e]['imprecise']==0:
					output+=['SVTYPE=INS;EVENT=%s;CHR2=%s;END=%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, sv_out[e]['gap'])]
				else:
					output+=['SVTYPE=INS;EVENT=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, \
					sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],sv_out[e]['CIpos2'][1], sv_out[e]['gap'])]
			elif(sv_out[e]['SVTYPE']=='DEL'):
				output+=['<DEL>','.','.']
				if sv_out[e]['imprecise']==0:
					output+=['SVTYPE=DEL;EVENT=%s;CHR2=%s;END=%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, abs(pos12.pos-pos11.pos)+1)]
				else:
					output+=['SVTYPE=DEL;EVENT=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, \
					sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],sv_out[e]['CIpos2'][1], abs(pos12.pos-pos11.pos)+1)]
			elif(sv_out[e]['SVTYPE']=='DUP'):
				output+=['<DUP>','.','.']
				if sv_out[e]['imprecise']==0:
					output+=['SVTYPE=DUP;EVENT=%s;CHR2=%s;END=%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, abs(pos12.pos-pos11.pos)+1)]
				else:
					output+=['SVTYPE=DUP;EVENT=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s;SVLEN=%s' % (event, pos12.chr, pos12.pos, \
					sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],sv_out[e]['CIpos2'][1], abs(pos12.pos-pos11.pos)+1)]
			else:
				SVTYPE_flag='BND'
				if len(sv_out.keys())==1:
					SVTYPE_flag='TRA'
				alt=''
				if sv_out[e]['STRAND'][1]=='+':
					alt=']'+sv_out[e]['pos2'].all+']'
				else:
					alt='['+sv_out[e]['pos2'].all+'['
				if sv_out[e]['STRAND'][0]=='+':
					if 'gap' in sv_out[e]:
						alt='N<INS>'+alt
					else:
						alt='N'+alt
				else:
					if 'gap' in sv_out[e]:
						alt=alt+'<INS>N'
					else:
						alt=alt+'N'
				output+=[alt,'.','.']
				if sv_out[e]['imprecise']==0:
					if 'gap' in sv_out[e]:
						output+=['SVTYPE=%s;EVENT=%s;STRANDS=%s;CHR2=%s;END=%s;INSLEN=%s' % (SVTYPE_flag,event,sv_out[e]['STRAND'],pos12.chr, \
						pos12.pos,sv_out[e]['gap'])]
					else:
						output+=['SVTYPE=%s;EVENT=%s;STRANDS=%s;CHR2=%s;END=%s' % (SVTYPE_flag,event,sv_out[e]['STRAND'],pos12.chr, pos12.pos)]
				else:
					if 'gap' in sv_out[e]:
							output+=['IMPRECISE;SVTYPE=%s;EVENT=%s;STRANDS=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s;INSLEN=%s' % (SVTYPE_flag,\
							event,sv_out[e]['STRAND'],pos12.chr, pos12.pos, sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],\
							sv_out[e]['CIpos2'][1],sv_out[e]['gap'])]
					else:
						output+=['IMPRECISE;SVTYPE=%s;EVENT=%s;STRANDS=%s;CHR2=%s;END=%s;CIPOS95=%s,%s;CIEND95=%s,%s' % (SVTYPE_flag,event,\
						sv_out[e]['STRAND'],pos12.chr, pos12.pos, sv_out[e]['CIpos1'][0], sv_out[e]['CIpos1'][1], sv_out[e]['CIpos2'][0],\
						sv_out[e]['CIpos2'][1])]
			output+=['GT:SU:SR','./.:%d:%d' % (sv_out[e]['SU'], sv_out[e]['SR'])]
			outf.write('\t'.join([str(x) for x in output])+'\n') 
			print_idx+=1
	outf.close()
