#!/usr/bin/python3
import sys
import os
import re
from math import log
from classes import *
from utils import function as f

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NaSV

read={}
sv_1={}
sv_2={}

def parse_maf():
	"""description"""
	read_maf(NaSV.opt_maf)
	convert_read()
	make_breakends()

def read_maf(maf):
	"""description"""
	inf=open(maf,'r')
	for l1 in inf:
		if l1[0]=='#' or l1=='\n':
			continue
		l2=next(inf)
		l3=next(inf)
		l4=next(inf)
		list1=l1.strip().split()
		list2=l2.strip().split()
		list3=l3.strip().split()
		list4=l4.strip().split()
		eg2=float(list1[2].split('=')[1])
		if -log(eg2)<NaSV.mq_cutoff:
			continue
		ref_chr, ref_st, ref_align_len, read_st=[list2[1], int(list2[2]), int(list2[3]), int(list3[2])]
		read_align_len, strand, read_len, read_id=[int(list3[3]), list3[4], int(list3[5]), '|'.join([list3[1],list3[5]])]
		read_ed=read_st+read_align_len-1
		ref_ed=ref_st+ref_align_len-1
		if strand=='-':
			read_ed=read_len-1-read_st
			read_st=read_ed-read_align_len+1
			ref_st,ref_ed=ref_ed,ref_st
		if read_id not in read:
			read[read_id]={}
		read[read_id][read_st]=READ(read_st,read_ed,ref_chr,ref_st,ref_ed,float(f.GC_ratio(list3[6])),eg2)
	inf.close()


def convert_read():
	"""description"""
	for i in read.keys():
		idx_list_raw=read[i].keys()
		if len(idx_list_raw)<2:#filter out non-fusion reads
			continue
		sorted(idx_list_raw,key = int)
		idx_list=[]
    #decide the first element and assemble from first element 
		len_i=int(i.split('|')[1])
		if list(idx_list_raw)[0] <= NaSV.start_cutoff*len_i:
			idx_list=list(idx_list_raw)
		else:
			continue
		key_=[]#key_=[breakend1_pos_gap_breakend2_pos,breakend3_pos_gap_breakend4_pos]
		value_=[]####value_=[[segment1_st,...],[segment1_ed,...][gap1,...][segment2_st,...],[segment2_ed,...]...]
		for j in range(0,len(idx_list)-1): 
			r1=read[i][idx_list[j]]
			r2=read[i][idx_list[j+1]]
			key=Breakend(ChrPos(r1.chr,r1.Red),r1.ori,r2.Qst-r1.Qed-1,ChrPos(r2.chr,r2.Rst),r2.ori)
			if len(idx_list)==2:
				value=[[ChrPos(r1.chr,r1.Rst)],[ChrPos(r1.chr,r1.Red)],[r2.Qst-r1.Qed-1],[ChrPos(r2.chr,r2.Rst)],[ChrPos(r2.chr,r2.Red)]]
				key_=[key]
				value_=value
				break
			elif j==0:
				value=[[ChrPos(r1.chr,r1.Rst)],[ChrPos(r1.chr,r1.Red)],[r2.Qst-r1.Qed-1],[ChrPos(r2.chr,r2.Rst)]]
			elif j==len(idx_list)-2:
				value=[[ChrPos(r1.chr,r1.Red)],[r2.Qst-r1.Qed-1],[ChrPos(r2.chr,r2.Rst)],[ChrPos(r2.chr,r2.Red)]]
			else:
				value=[[ChrPos(r1.chr,r1.Red)],[r2.Qst-r1.Qed-1],[ChrPos(r2.chr,r2.Rst)]]
			key_.append(key)
			value_=value_+value
		if (key_[0].ori1==key_[-1].ori2=='-' or (f.key_symmetry(key_) and value_[0][0].pos>value_[-1][0].pos)):
			key_ = f.key_reverse(key_)
			value_ = list(reversed(value_))
		if len(sv_1.keys())==0:
			sv_1['|'.join([x.all for x in key_])]=value_
		else:
			merge_=0
			for k in sv_1.keys():
				k_list=k.split('|')
                ###if intersect,merge
				if(len(key_)==len(k_list) and f.breakends_intersect([x.all for x in key_],k_list)):
					for l in range(0,len(sv_1[k])):
						sv_1[k][l].append(value_[l][0])
					merge_=1
					break
			if merge_==0:
				sv_1['|'.join([x.all for x in key_])]=value_


def make_breakends():
	"""description"""
	for i in sv_1.keys():
		if len(sv_1[i][0])<NaSV.read_cutoff:
			continue
		pos_list=[]
        #final result: 1)mode if count(mode)/len(array)>=50%;2)count mean position 
		for j in sv_1[i]:
			pos,_=f.representative(j)
			pos_list.append(pos)
		final_key=''
        #get key_id elements
		count=0
		for j in i.split('|'):
			m=re.match(r'(\w+):\d+([-+])_-?\d+_(\w+):\d+([-+])',j)
			l_chr,l_ori,r_chr,r_ori=m.groups()
			if count==0:
				final_key = '%s%s_%s_%s%s' % (pos_list[3*count+1].all,l_ori,pos_list[3*count+2],pos_list[3*count+3].all,r_ori)
			else:
				final_key = '|'.join([final_key,'%s%s_%s_%s%s' % (pos_list[3*count+1].all,l_ori,pos_list[3*count+2],pos_list[3*count+3].all,r_ori)])
			count += 1
		sv_2[final_key]=[len(sv_1[i][1]),pos_list[0],pos_list[-1],i]