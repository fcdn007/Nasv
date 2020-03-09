#!/usr/bin/python3
import sys
import os
import re
from classes import *
from utils import function as f
from utils.parse_maf import sv_2
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NaSV

key_region={}
sv_3={}

def parse_breakpoints():
    make_keyRegion()
    make_sv()

def make_keyRegion():
	for i in sv_2.keys():
		list_=[int(x) if re.match('^-?\d+$',x) else ChrPos_plus(x) for x in re.split('\||_',i)]
		next_=0
		key_region[i]=[]
		key_region[i].append(ChrRegion(sv_2[i][1],list_[0]))
		for j in range(1,len(list_)-1):
			if next_:
				next_=0
				continue
			elif (j+1) % 3 == 0:
				next_=1
				key_region[i].append(ChrRegion(list_[j],list_[j+1]))
			elif list_[j]>NaSV.merge_cutoff:
				key_region[i].append(list_[j])
		key_region[i].append(ChrRegion(list_[-1],sv_2[i][2]))

def make_sv():
	for i in sorted(key_region.keys(),key=lambda x: len(x.split('|')), reverse=True):
		if i not in key_region:
			continue
		readNum=sv_2[i][0]
		sv_3[i]={}
		sv_3[i]['SWR']=readNum
		sv_3[i]['region']={}
		for x in range(len(key_region[i])):
			sv_3[i]['region'][x]={}
			if isinstance(key_region[i][x],int):
				sv_3[i]['region'][x][key_region[i][x]]=readNum
			else:
				sv_3[i]['region'][x][key_region[i][x].all]=readNum
		for j in sorted(key_region.keys(),key=lambda x: len(x.split('|')), reverse=True):
			if i == j or len(key_region[j])>len(key_region[i]):
				continue
			match=[]
			match_key=[]
			merge_, match_=[0, 0]
			for i_ in range(0,len(key_region[i])):
				i_end,j_end=[0,-1]
				if i_==0:
					i_end=-1
				elif i_==len(key_region[i])-1:
					i_end=1
				i_list=key_region[i][i_]
				j_list=key_region[j][0]
				match_=f.compare_region(j_list,i_list,j_end,i_end)
				if (match_==1 and i_<=len(key_region[i])-len(key_region[j])) or (match_==2 and i_>=len(key_region[j])-1):
					match_key.append(i_)
					match.append(match_)
			match_break=0
			for k in range(0,len(match)):
				if match_break:
					break
				match_=match[k]
				match_key_=match_key[k]
				if match_==1:
					for j_ in range(1,len(key_region[j])):
						i_end,j_end=[0,0]
						if j_==len(key_region[j])-1:
							j_end=1
						if match_key_+j_==len(key_region[i])-1:
							i_end=1
						match_=f.compare_region(key_region[j][j_],key_region[i][match_key_+j_],j_end,i_end)
						if match_==0 or match_==2:
							break
					if match_==1:
						readNum=sv_2[j][0]
						tmp='%s:NONE_%s%s' % (key_region[j][0].chr, key_region[j][0].pos2, key_region[j][0].ori)
						if tmp in sv_3[i]['region'][match_key_]:
							sv_3[i]['region'][match_key_][tmp]+=readNum
						else:
							sv_3[i]['region'][match_key_][tmp]=readNum
						for h in range(1,len(key_region[j])-1):
							if isinstance(key_region[j][h],int):
								if key_region[j][h] in sv_3[i]['region'][h+match_key_]:
									sv_3[i]['region'][h+match_key_][key_region[j][h]]+=readNum
								else:
									sv_3[i]['region'][h+match_key_][key_region[j][h]]=readNum
							else:
								if key_region[j][h].all in sv_3[i]['region'][h+match_key_]:
									sv_3[i]['region'][h+match_key_][key_region[j][h].all]+=readNum
								else:
									sv_3[i]['region'][h+match_key_][key_region[j][h].all]=readNum
						tmp='%s:%s_NONE%s' % (key_region[j][-1].chr, key_region[j][-1].pos1, key_region[j][-1].ori)
						if tmp in sv_3[i]['region'][len(key_region[j])-1+match_key_]:
							sv_3[i]['region'][len(key_region[j])-1+match_key_][tmp]+=readNum
						else:
							sv_3[i]['region'][len(key_region[j])-1+match_key_][tmp]=readNum
						del key_region[j]
						match_break=1
				else:
					for j_ in range(1,len(key_region[j])):
						i_end,j_end=[0,0]
						if j_==len(key_region[j])-1:
							j_end=1
						if match_key_-j_==0:
							i_end=-1
						match_=f.compare_region(key_region[j][j_],key_region[i][match_key_-j_],j_end,i_end)
						if match_==0 or match_==1:
							break
					if match_==2:
						readNum=sv_2[j][0]
						if key_region[j][0].ori=='+':
							t_ori='-'
						else:
							t_ori='+'
						tmp='%s:%s_NONE%s' % (key_region[j][0].chr, key_region[j][0].pos2, t_ori)
						if tmp in sv_3[i]['region'][match_key_]:
							sv_3[i]['region'][match_key_][tmp]+=readNum
						else:
							sv_3[i]['region'][match_key_][tmp]=readNum
						for h in range(1,len(key_region[j])-1):
							if f.regionKey_reverse(key_region[j][h]) in sv_3[i]['region'][match_key_-h]:
								sv_3[i]['region'][match_key_-h][f.regionKey_reverse(key_region[j][h])]+=readNum
							else:
								sv_3[i]['region'][match_key_-h][f.regionKey_reverse(key_region[j][h])]=readNum
						if key_region[j][-1].ori=='+':
							t_ori='-'
						else:
							t_ori='+'
						tmp='%s:NONE_%s%s' % (key_region[j][-1].chr, key_region[j][-1].pos1, t_ori)
						if tmp in sv_3[i]['region'][1-len(key_region[j])+match_key_]:
							sv_3[i]['region'][1-len(key_region[j])+match_key_][tmp]+=readNum
						else:
							sv_3[i]['region'][1-len(key_region[j])+match_key_][tmp]=readNum
						del key_region[j]
						match_break=1
