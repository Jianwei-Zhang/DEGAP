import re
import os
import sys
import getopt
import pysam
from pysam import AlignmentFile
import math
import Bio
from Bio import SeqIO
import GapFiller
from GapFiller import GapFiller
import CtgLinker
from CtgLinker import CtgLinker
import selectRawReads
from selectRawReads import selectRawReads

def usage():
	print ("--reads HiFi_reads.fasta")
	print ("-o | --out ./path/")
	print ("-t | --thread thread number")
	print ("--remove 1 | 2 | 3 default:2 1:only keep final result; 2: keep every round basic result ; 3 : keep all files")
	print ("--edge Edge Controller set max Edge length(missequening)")
	print ("--filterDepth num default:None. You can filtered HiFi reads by mapped depth on contig set. if num==0.3 means: mapped Hifi reads on depth>=0.3*avgdepth and depth<=(2-0.2)*avgdepth will be filtered and will not be used in whole project")
	print ("--mode gapfiller | ctglinker")
	print ("\n\ngapfiller\n")
	print ("\t--seqleft sequence before GAP")
	print ("\t--seqright sequence after GAP")
	print ("\t--flag left | right choose seqleft or seqright as seed to fill the gap(default s)")
	print ("\n\nctglinker\n")
	print ("\t--ctgseq contig set")
	print ("-h | --help help")

def getoptions():
	try:
		opts,args= getopt.getopt(sys.argv[1:], "o:t:h",["reads=","seqleft=","seqright=","flag=","out=","ctgseq=","mode=","edge=","thread=","help=","remove=","filterDepth="])
	except getopt.GetoptError:
		print ("Apeared Error Parameter!!")
		usage()
		sys.exit()
	
	if len(opts)==0:
		print ("NO PARAMETER!!!")
		print ("Use '-h | --help' for same information")
		sys.exit()
	
	mode=''
	reads=''
	out=''
	edge=500
	filterDepth=None
	thread='20'
	remove=2
	
	for opt,value in opts:
		if opt in ("--mode"):
			mode=value
		elif opt in ("-t","--thread"):
			thread=str(int(value))
		elif opt in ("--reads"):
			reads=value
		elif opt in ("-o","--out"):
			out=value
		elif opt in ("-h","--help"):
			usage()
			sys.exit()
		elif opt in ("--edge"):
			edge=int(value)
		elif opt in ("--filterDepth"):
			filterDepth=float(value)
		elif opt in ("--remove"):
			remove=int(value)
	if mode=="gapfiller":
		seqleft=""
		seqright=""
		flag='left'
		for opt,value in opts:
			if opt in ("--seqleft"):
				if os.path.exists(value)==True and os.path.getsize(value)!=0:
					seqleft=value
				else:
					print ("seqleft file don't exist or empty!")
					sys.exit()
			elif opt in ("--seqright"):
				if os.path.exists(value)==True and os.path.getsize(value)!=0:
					seqright=value
				else:
					print ("seqright file  don't exist or empty!")
					sys.exit()
			elif opt in ("--flag"):
				flag=value
		return [mode,remove,thread,reads,out,seqleft,seqright,flag,edge,filterDepth]
	elif mode=="ctglinker":
		seqfile=""
		for opt,value in opts:
			if opt in ("--ctgseq"):
				if os.path.exists(value)==True and os.path.getsize(value)!=0:
					seqfile=value
				else:
					print ("contigs file don't exist or empty")
					sys.exit()
		return [mode,remove,thread,reads,out,seqfile,edge,filterDepth]
	else:				
		print ("You should use gapfiller or ctglinker!")
		sys.exit()
parameter=getoptions()

reads=parameter[3]
out=parameter[4]

#make outpit file
if out[-1]=="/":
	out=out[:-1]
	parameter[4]=out
if os.path.exists(out)!=True:
	os.makedirs(out)

#make reads index
print ("BUILD RAWREADS DICT")
readsdict=SeqIO.index_db(out+"/reads.idx",reads,"fasta")
print ("BUILD DICT SUCCEED")
parameter.append(readsdict)

pwd1=out+"/HiFi.reads.stat"
if os.path.exists(pwd1)==True and os.path.getsize(pwd1)!=0:
	file1=open(pwd1,'r')
	for i in file1:
		i1=i.rstrip().split("\t")
		if i1[0]=="MaxLength":
			lenmax=int(i1[1])
		elif i1[0]=="SeedLength":
			seedlen=int(float(i1[1]))
else:
	lenmax=0
	file1=open(pwd1,'w')
	n=0
	total_len=0
	for gseq in SeqIO.parse(reads,'fasta'):
		n=n+1
		if len(gseq.seq)>lenmax:
			lenmax=len(gseq.seq)
		total_len=total_len+len(gseq.seq)
	l='Number\t'+str(n)+"\nTolalLenth\t"+str(total_len)+"\nMaxLength\t"+str(lenmax)+"\n"
	file1.writelines(l)
	a=10**(int(math.log(lenmax,10)))
	b=lenmax/a+1
	seedlen=a*b
	l='SeedLength\t'+str(seedlen)+"\n"
	file1.writelines(l)

file1.close()
print (lenmax)
print (seedlen)
if parameter[-2]!=None:
	selectedReads=selectRawReads(parameter,seedlen)

	parameter=parameter[:-1]

	print ("BUILD RAWREADS DICT")
	readsdict2=SeqIO.index_db(out+"/selectedReads.idx",selectedReads.readFile,"fasta")
	print ("BUILD DICT SUCCEED")
	parameter.append(readsdict2)
	parameter[3]=selectedReads.readFile

parameter.append(lenmax)
parameter.append(seedlen)
print (parameter)
if parameter[0]=="gapfiller":
	GapFiller(parameter)
else:
	CtgLinker(parameter)
print ('welldone')


