#!/usr/bin/python3
import sys
import os
import re
from math import log,ceil,sqrt
from scipy import stats
import numpy as np
from classes import ChrPos, ChrPos_plus, Breakend, Breakend_plus
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
import NaSV

def GC_ratio(seq):
    len_=0
    num=0
    for i in seq[:]:
        if i=='-':
            continue
        else:
            len_+=1
        if re.match(r'C|G|c|g',i):
            num+=1
    return '%.4f' % (float(num)/len_)

def breakends_intersect(l1,l2):
    merge_=0
    for i in range(0,len(l1)):
        m1=Breakend_plus(l1[i])
        m2=Breakend_plus(l2[i])
        if m1.chr1 == m2.chr1 and m1.ori1== m2.ori1 and m1.chr2 == m2.chr2 and m1.ori2== m2.ori2 and abs(m1.pos1-m2.pos1)<=NaSV.merge_bp_cutoff \
           and abs(m1.pos2-m2.pos2)<=NaSV.merge_bp_cutoff and abs(m1.gap-m2.gap)<=2*NaSV.merge_bp_cutoff:
            merge_=1
        else:
            merge_=0
            break
    return merge_

def reverse_ori(l1):
    v=[]
    for ori in l1:
        if ori == '-':
            v.append('+')
        else:
            v.append('-')
    return v

def key_reverse(l1):
    l2=[]
    for i in reversed(l1):     
        sign_r,sign_l=reverse_ori([i.ori1,i.ori2]) 
        l2.append(Breakend(i.e2,sign_l,i.gap,i.e1,sign_r))
    return l2

def which_tpye(e):
    if isinstance(e, int):
        return 2 #gap
    else:
        return 1 #region_coordiante

def compare_region(e1,e2,end1,end2):
#e1:small region; e2:large region
#0:not overlap; 1:positive overlap; 2: negative overlap; 3: gap overlap
    e1_type=which_tpye(e1)
    e2_type=which_tpye(e2)
    if end2==0:
        if e1_type==e2_type:
            if e1_type==1:
                if e1.ori==e2.ori=="+" and e1.chr == e2.chr and e1.pos1>=e2.pos1-NaSV.merge_cutoff and e1.pos2<=e2.pos2+NaSV.merge_cutoff:
                    return 1
                elif e1.ori==e2.ori=="-" and e1.chr == e2.chr and e1.pos1<=e2.pos1+NaSV.merge_cutoff and e1.pos2>=e2.pos2-NaSV.merge_cutoff:
                    return 1
                elif e1.ori=="-" and e2.ori=="+" and e1.chr == e2.chr and e1.pos2>=e2.pos1-NaSV.merge_cutoff and e1.pos1<=e2.pos2+NaSV.merge_cutoff:
                    return 2
                elif e1.ori=="+" and e2.ori=="-" and e1.chr == e2.chr and e1.pos2<=e2.pos1+NaSV.merge_cutoff and e1.pos1>=e2.pos2-NaSV.merge_cutoff:
                    return 2
                else:
                    return 0
            else:
                if int(e1)<=int(e2):
                    return 3
                else:
                    return 0
        elif (e1_type==2 and e2_type==1):
            if abs(e2.pos1-e2.pos2)+1>=e1-NaSV.merge_cutoff:
                return 3
        elif (e1_type==1 and e2_type==2):
            if abs(e1.pos1-e1.pos2)+1<=e2+NaSV.merge_cutoff:
                return 3
        elif (e1_type==e2_type==2):
            if abs(e1-e2)+1<=NaSV.merge_cutoff:
                return 3
        else:
            return 0
    elif end1==0 and end2!=0:
        return 0
    elif end1==-1 and end2==-1:
        if e2.ori==e1.ori and e2.chr==e1.chr and abs(e2.pos2-e1.pos2)<=NaSV.merge_cutoff: 
            return 1       
    elif end1==-1 and end2==1:
        if e2.ori!=e1.ori and e2.chr==e1.chr and abs(e2.pos1-e1.pos2)<=NaSV.merge_cutoff: 
            return 2        
    elif end1==1 and end2==-1:
        if e2.ori!=e1.ori and e2.chr==e1.chr and abs(e2.pos2-e1.pos1)<=NaSV.merge_cutoff: 
            return 2         
    elif end1==1 and end2==1:
        if e2.ori==e1.ori and e2.chr==e1.chr and abs(e2.pos1-e1.pos1)<=NaSV.merge_cutoff: 
            return 1
        
def key_symmetry(l):
    l_t=[]
    for i in l:
        l_t=l_t+i.all.split('_')
    result=0
    for i in range(int(len(l_t)/2)):
        if i %3==1:
            if abs(int(l_t[i])-int(l_t[-int(i+1)]))<=NaSV.merge_cutoff:
                result=1
            else:
                result=0
        else:
            e1=l_t[i].split(":")
            e2=l_t[-int(i+1)].split(":")
            if e1[0]==e2[0] and e1[1][-1]!=e2[1][-1] and abs(int(e1[1][:-1])-int(e2[1][:-1]))<=NaSV.merge_cutoff:
                result=1
            else:
                reslut=0
        if result==0:
            break
    return result

def CI(l1,ratio=0.95):
	l2=l1
	if not isinstance(l1[0],int):
		l2=[x.pos for x in l1]
	L=np.array(l2)
	L_mean=np.mean(L)
	L_SE=np.std(L)/np.sqrt(L.size)
	if L_SE==0:
		return [0,0]
	else:
		up,down=stats.norm.interval(ratio,loc=L_mean,scale=L_SE)-L_mean
		return ["%.2f" % up.item(), "%.2f" % down.item()]
    
def representative(l):
    sum_=0
    count={}
    precise=0
    for k in l:
        value=''
        if isinstance(k,int):
            value=int(k)
        else:
            value=k.pos
        sum_ += value
        if value in count:
            count[value]+=1
        else:
            count[value]=1
    mode,mode_ratio=sorted(count.items(),key=lambda x:x[1])[-1]
    if mode_ratio/float(len(l))>=NaSV.precise_cutoff:
        pos = int(mode)
        precise=1       
    else:
        pos = int(sum_/len(l))
    if not isinstance(l[0],int):
        pos = ChrPos(l[0].chr,pos)
    return pos,precise

def regionKey_reverse(s):
    if isinstance(s,int):
        return s
    else:
        l=s.all.split('_')
        ori='+'
        l_chr, l_pos1=l[0].split(':')
        l_pos2=l[1][:-1]
        if l[1][-1]=='+':
            ori='-'
        return '%s:%s_%s%s' % (l_chr,l_pos2,l_pos1,ori)
