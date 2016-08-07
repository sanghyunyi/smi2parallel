#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@brief convert .smi subtitle file into aligned parallel corpora text file using dp 
license: GPL

@version: 1.0.0
@author: SangHyun Yi <sangyi0914@gmail.com>
'''
__author__ = "SangHyun Yi <sangyi0914@gmail.com>"
__date__ = "2016/06/15"
__version__ = "1.0.0"
__version_info__ = (1, 0, 0)

###################################################################################################
import chardet #@UnresolvedImport
import re
import sys
import codecs
import os
import operator
import io
import math
import numpy as np

###################################################################################################
def usage(msg=None, exit_code=1):
	print_msg = """
usage %s smifile.smi [...]
	convert smi into srt subtitle file with same filename.
	By MoonChang Chae <mcchae@gmail.com>
""" % os.path.basename(sys.argv[0])
	if msg:
		print_msg += '%s\n' % msg
	print print_msg
	sys.exit(exit_code)

class Alignment(object):
    def __init__(self, x, y, d):
        self.x = x
        self.y = y
        self.d = d

def getOverlap(a, b):
    return(min(a[1], b[1]) - max(a[0], b[0]))

def TwoSideDistance(e, k, en, kn): #get list of e and k sentences
    constant = 10000
    if en == 0 or kn == 0:
        return(0)
    else:
        overlap = getOverlap((e[0].start_ms, e[en-1].end_ms),(k[0].start_ms, k[kn-1].end_ms))-constant*abs(en-kn)
        return(overlap)
    
def SeqAlign(x, y): #x,y are list of sentences to be aligned. nx, ny are number of sentences to be aligned
    nx = len(x)
    ny = len(y)
    ratio = float(nx)/float(ny)
    distances = np.zeros((nx + 1, ny + 1))
    path_x = np.zeros((nx + 1, ny + 1))
    path_y = np.zeros((nx + 1, ny + 1))
    for j in range(ny + 1):
        ipivot = int(j*ratio)
        istart = max(0,ipivot-30)
        iend = min(nx+1,ipivot+30)
        for i in range(istart,iend):
            d1 = distances[i-1, j-1] + TwoSideDistance([x[i-1]],[y[j-1]],1,1) if i > 0 and j > 0 else -sys.maxint
            d2 = distances[i-1, j] + TwoSideDistance([x[i-1]],[y[j-1]],1,0) if i > 0 else -sys.maxint
            d3 = distances[i, j-1] + TwoSideDistance([x[i-1]],[y[j-1]],0,1) if j > 0 else -sys.maxint

            d4 = distances[i-2, j-1] + TwoSideDistance([x[i-2],x[i-1]],[y[j-1]],2,1) if i > 1 and j > 0 else -sys.maxint
            d5 = distances[i-1, j-2] + TwoSideDistance([x[i-1]],[y[j-2],y[j-1]],1,2) if i > 0 and j > 1 else -sys.maxint
            d6 = distances[i-2, j-2] + TwoSideDistance([x[i-2],x[i-1]],[y[j-2],y[j-1]],2,2) if i > 1 and j > 1 else -sys.maxint

            d7 = distances[i-3, j-1] + TwoSideDistance([x[i-3],x[i-2],x[i-1]],[y[j-1]],3,1) if i > 2 and j > 0 else -sys.maxint
            d8 = distances[i-3, j-2] + TwoSideDistance([x[i-3],x[i-2],x[i-1]],[y[j-2],y[j-1]],3,2) if i > 2 and j > 1 else -sys.maxint
            d9 = distances[i-1, j-3] + TwoSideDistance([x[i-1]],[y[j-3],y[j-2],y[j-1]],1,3) if i > 0 and j > 2 else -sys.maxint
            d10 = distances[i-2, j-3] + TwoSideDistance([x[i-2],x[i-1]],[y[j-3],y[j-2],y[j-1]],2,3) if i > 1 and j > 2 else -sys.maxint
            d11 = distances[i-3, j-3] + TwoSideDistance([x[i-3],x[i-2],x[i-1]],[y[j-3],y[j-2],y[j-1]],3,3) if i > 2 and j > 2 else -sys.maxint

            dmax = max([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11])
            if dmax == -sys.maxint:
                distances[i, j] = 0
            elif dmax == d1:
                distances[i, j] = d1
                path_x[i, j] = i - 1
                path_y[i, j] = j - 1
            elif dmax == d2:
                distances[i, j] = d2
                path_x[i, j] = i - 1
                path_y[i, j] = j
            elif dmax == d3:
                distances[i, j] = d3
                path_x[i, j] = i
                path_y[i, j] = j - 1
            elif dmax == d4:
                distances[i, j] = d4
                path_x[i, j] = i - 2
                path_y[i, j] = j - 1
            elif dmax == d5:
                distances[i, j] = d5
                path_x[i, j] = i - 1
                path_y[i, j] = j - 2
            elif dmax == d6:
                distances[i, j] = d6
                path_x[i, j] = i - 2
                path_y[i, j] = j - 2
            elif dmax == d7:
                distances[i, j] = d7
                path_x[i, j] = i - 3
                path_y[i, j] = j - 1
            elif dmax == d8:
                distances[i, j] = d8
                path_x[i, j] = i - 3
                path_y[i, j] = j - 2
            elif dmax == d9:
                distances[i, j] = d9
                path_x[i, j] = i - 1
                path_y[i, j] = j - 3
            elif dmax == d10:
                distances[i, j] = d10
                path_x[i, j] = i - 2
                path_y[i, j] = j - 3
            else:
                distances[i, j] = d11
                path_x[i, j] = i - 3
                path_y[i, j] = j - 3

    Alignlist = []
    i = nx
    j = ny
    while i > 0 or j > 0:
        oi = path_x[i, j]
        oj = path_y[i, j]
        di = i - oi
        dj = j - oj
        if di == 1 and dj == 1:
            ralign = Alignment([x[i-1]], [y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 1 and dj == 0:
            ralign = Alignment([x[i-1]], [], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 0 and dj == 1:
            ralign = Alignment([], [y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 2 and dj == 1:
            ralign = Alignment([x[i-2], x[i-1]], [y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 1 and dj == 2:
            ralign = Alignment([x[i-1]], [y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 2 and dj == 2:
            ralign = Alignment([x[i-2], x[i-1]], [y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)

        elif di == 1 and dj == 3:
            ralign = Alignment([x[i-1]], [y[j-3], y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 2 and dj == 3:
            ralign = Alignment([x[i-2], x[i-1]], [y[j-3], y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 3 and dj == 1:
            ralign = Alignment([x[i-3], x[i-2], x[i-1]], [y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        elif di == 3 and dj == 2:
            ralign = Alignment([x[i-3], x[i-2], x[i-1]], [y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)
        else:
            ralign = Alignment([x[i-3], x[i-2], x[i-1]], [y[j-3], y[j-2], y[j-1]], distances[i, j] - distances[i-di, j-dj])
            Alignlist.append(ralign)

        i = int(oi)
        j = int(oj)
    Alignlist.reverse()
    return(Alignlist)

###################################################################################################
class smiItem(object):
	def __init__(self):
		self.start_ms = 0L
		self.start_ts = '00:00:00,000'
		self.end_ms = 0L
		self.end_ts = '00:00:00,000'
		self.contents = None
		self.linecount = 0
                self.isko = 1 # ko = 1, en = 0
	@staticmethod
	def ms2ts(ms):
		hours = ms / 3600000L
		ms -= hours * 3600000L
		minutes = ms / 60000L
		ms -= minutes * 60000L
		seconds = ms / 1000L
		ms -= seconds * 1000L
		s = '%02d:%02d:%02d,%03d' % (hours, minutes, seconds, ms)
		return s
	def convertSrt(self):
		if self.linecount == 4:
			i=1 #@UnusedVariable
		# 2) remove new-line
		self.contents = re.sub(r'\s+', ' ', self.contents)
		# 3) remove web string like "&nbsp";
		self.contents = re.sub(r'&[a-z]{2,5};', '', self.contents)
		# 4) replace "<br>" with ' ';
		self.contents = re.sub(r'(<br>)+', ' ', self.contents, flags=re.IGNORECASE)
                self.contents = re.sub(r'(<i>)+', '', self.contents, flags=re.IGNORECASE)
                self.contents = re.sub(r'(</i>)+', '', self.contents, flags=re.IGNORECASE)
		# 5) find all tags
		fndx = self.contents.find('<')
		if fndx >= 0:
			contents = self.contents
			sb = self.contents[0:fndx]
			contents = contents[fndx:]
			while True:
				m = re.match(r'</?([a-z]+)[^>]*>([^<>]*)', contents, flags=re.IGNORECASE)
				if m == None: break
				contents = contents[m.end(2):]
				#if m.group(1).lower() in ['font', 'b', 'i', 'u']:
				if m.group(1).lower() in ['b', 'i', 'u']:
					sb += m.string[0:m.start(2)]
				sb += m.group(2)
			self.contents = sb
		self.contents = self.contents.strip()
		self.contents = self.contents.strip('\n')
	def __repr__(self):
		s = '%d:%d:<%s>:%d' % (self.start_ms, self.end_ms, self.contents, self.linecount)
		return s

###################################################################################################
def convertSMI(smi_file):
	if not os.path.exists(smi_file):
		sys.stderr.write('Cannot find smi file <%s>\n' % smi_file)
		return False
	rndx = smi_file.rfind('.')
	srt_file = '%s.txt' % smi_file[0:rndx]

	ifp = open(smi_file)
	smi_sgml = ifp.read()#.upper()
	ifp.close()
	chdt = chardet.detect(smi_sgml)
	print(chdt['encoding'])
	if chdt['encoding'] == 'windows-1252':
            print("ge")
            smi_sgml = ''
        elif chdt['encoding'] == 'ISO-8859-2':
            print("ck")
            smi_sgml = ''
        elif chdt['encoding'] == None:
            print("dg")
            smi_sgml = ''
        else:
            smi_sgml = unicode(smi_sgml, chdt['encoding']).encode('utf-8')            

	# skip to first starting tag (skip first 0xff 0xfe ...)
	try:
		fndx = smi_sgml.find('<SYNC')
	except Exception, e:
		print chdt
		raise e
	if fndx < 0:
		return False
	smi_sgml = smi_sgml[fndx:]
	lines = smi_sgml.split('\n')
        ko_list = []
        en_list = []
	sync_cont = ''
	si = None
	last_si = None
	linecnt = 0
	for line in lines:
		linecnt += 1
		sndx = line.upper().find('<SYNC')
		if sndx >= 0:
			m = re.search(r'<sync\s+s\s*t\s*a\s*r\s*t\s*=\s*(-?\d+)>\s*<\s*P\s*C\s*l\s*a\s*s\s*s\s*=\s*(\w+)\s*>(.*)$', line, flags=re.IGNORECASE)
			if not m:
				raise Exception('Invalid format tag of <Sync start=nnnn> with "%s"' % line)
			sync_cont += line[0:sndx]
			last_si = si
			if last_si != None:
				last_si.end_ms = long(m.group(1))
                                if re.search(r'K\s*R\s*C\s*C\s*', m.group(2).upper()):
                                    last_si.isko = 1
                                elif re.search(r'S\s*U\s*B\s*T\s*T\s*L', m.group(2).upper()):
                                    last_si.isko = 1
                                else:
                                    last_si.isko = 0
				last_si.contents = sync_cont
				last_si.linecount = linecnt
                                if last_si.isko:
                                    ko_list.append(last_si)
                                else:
                                    en_list.append(last_si)
				#print '[%06d] %s' % (linecnt, last_si)
			sync_cont = m.group(3)
			si = smiItem()
			si.start_ms = long(m.group(1))
		else:
			sync_cont += line
        
        ko_list.sort(key = lambda x: x.start_ms )
        en_list.sort(key = lambda x: x.start_ms )
        for si in ko_list:
            si.convertSrt()
        for si in en_list:
            si.convertSrt()
        print "Alignment Start!"
        aligned = SeqAlign(en_list,ko_list)
        print "Alignment Complete!"
        ofp = open(srt_file, 'w')

        for alignment in aligned:
            xtext = ''
            ytext = ''
            for x in alignment.x:
                if len(x.contents) > 0:
                    xtext += x.contents + ' '
            for y in alignment.y:
                if len(y.contents) > 0:
                    ytext += y.contents + ' '
            if len(xtext) > 0 and len(ytext) > 0:
                outtext = xtext + '\n' + ytext + '\n'
                ofp.write(outtext)
	ofp.close()
	return True

###################################################################################################
def doConvert():
	if len(sys.argv) <= 1:
		usage()
	for smi_file in sys.argv[1:]:
                print(smi_file)
		if convertSMI(smi_file):
			print "Converting <%s> OK!" % smi_file
		else:
			print "Converting <%s> Failure!" % smi_file
			
	
###################################################################################################
if __name__ == '__main__':
	doConvert()
