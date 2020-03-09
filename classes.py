#!/usr/bin/python3
import re

class ChrPos:
    def __init__(self, chr_, pos_, ori='N'):
        self.chr=chr_
        self.pos=int(pos_)
        self.ori=ori
        if self.ori=='N':
            self.all=':'.join([chr_,str(pos_)])
        else:
            self.all='%s:%s%s' % (chr_,pos_,ori)

class ChrPos_plus:
	def __init__(self, s):
		m=re.match('(chr.+):(\d+)',s)
		self.chr,pos_=m.groups()
		self.pos=int(pos_)
		self.all=s
		if s[-1]=='+' or s[-1]=='-':
			self.ori=s[-1]

class ChrRegion: 
    def __init__(self, ChrPos1, ChrPos2):
        if isinstance(ChrPos1,str):
            chr_,pos_=ChrPos1.split(":")
            ChrPos1=ChrPos(chr_,pos_)
        if isinstance(ChrPos2,str):
            chr_,pos_=ChrPos2.split(":")
            ChrPos2=ChrPos(chr_,pos_)
        self.chr=ChrPos1.chr
        self.pos1=ChrPos1.pos
        self.pos2=ChrPos2.pos
        if self.pos1>self.pos2:
            self.ori='-'
        elif self.pos1<self.pos2:
            self.ori='+'
        self.all='%s:%s_%s%s' % (ChrPos1.chr,ChrPos1.pos,ChrPos2.pos,self.ori)
    
class Breakend:
    def __init__(self, ChrPos1, strand1, num, ChrPos2, strand2):
        if isinstance(ChrPos1,str):
            chr_,pos_=ChrPos1.split(":")
            ChrPos1=ChrPos(chr_,pos_)
        if isinstance(ChrPos2,str):
            chr_,pos_=ChrPos2.split(":")
            ChrPos2=ChrPos(chr_,pos_)
        self.e1=ChrPos1
        self.e2=ChrPos2
        self.chr1=ChrPos1.chr
        self.chr2=ChrPos2.chr
        self.pos1=ChrPos1.pos
        self.pos2=ChrPos2.pos
        self.gap=num
        self.ori1=strand1
        self.ori2=strand2
        self.all='%s%s_%s_%s%s' % (ChrPos1.all,strand1,num,ChrPos2.all,strand2)

class Breakend_plus:
	def __init__(self, s):
		m=re.match(r'(\w+):(\d+)([-+])_(-?\d+)_(\w+):(\d+)([-+])',s)
		chr1,pos1,ori1,gap,chr2,pos2,ori2=m.groups()	
		self.e1=ChrPos(chr1,pos1)
		self.e2=ChrPos(chr2,pos2)
		self.chr1=chr1
		self.chr2=chr2
		self.pos1=int(pos1)
		self.pos2=int(pos2)
		self.gap=int(gap)
		self.ori1=ori1
		self.ori2=ori2
		self.all=s

class READ:
	def __init__(self, read_st,read_ed,ref_chr,ref_st,ref_ed,gc,eg2):
		self.Qst=read_st
		self.Qed=read_ed
		self.chr=ref_chr
		self.Rst=ref_st
		self.Red=ref_ed
		self.gc=gc
		self.e=eg2
		if self.Red>self.Rst:
			self.ori='+'
		elif self.Rst>self.Red:
			self.ori='-'
		self.all="%s_%s|%s:%s-%s|%s|%s" % (read_st,read_ed,ref_chr,ref_st,ref_ed,gc,eg2)
