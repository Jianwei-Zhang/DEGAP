import re
import sys
import os
import operator
import numpy as np
import Bio
from Bio import SeqIO
import pysam
from pysam import AlignmentFile

class selectRawReads(object):
	def __init__(self,parameterlist,seedLen):
		self.seedLen=seedLen
		if len(parameterlist)==9:
			self.mode,self.remove,self.thread,self.reads,self.out,self.genomeSeq,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict=parameterlist
			self.inputSeq=self.genomeSeq
			self.inputSeqDict={}
			for gseq in SeqIO.parse(self.inputSeq,'fasta'):
				self.inputSeqDict[gseq.id]=len(gseq.seq)
		else:
			self.mode,self.remove,self.thread,self.reads,self.out,self.seqLeft,self.seqRight,self.flag,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict=parameterlist
			self.inputSeq=self.out+"/Genome.inputCtg.fa"
			filet=open(self.inputSeq,'w')
			self.inputSeqDict={}
			for gseq in SeqIO.parse(self.seqLeft,'fasta'):
				l='>'+gseq.id+"\n"+gseq.seq+"\n"
				filet.writelines(l)
				self.inputSeqDict[gseq.id]=len(gseq.seq)
			for gseq in SeqIO.parse(self.seqRight,'fasta'):
				l='>'+gseq.id+"\n"+gseq.seq+"\n"
				filet.writelines(l)
				self.inputSeqDict[gseq.id]=len(gseq.seq)
			filet.close()
		self.selectRawReadsLog=self.out+"/Genome.inputCtg.mappedReads.log"
		self.alnTag=1	
		if os.path.exists(self.selectRawReadsLog)==True and os.path.getsize(self.selectRawReadsLog)!=0:
			file1=open(self.selectRawReadsLog,'r')
			for row in file1:
				row1=row.rstrip().split('\t')
				if row1[0]=='alnFile':
					self.alnFile=row1[1]
				elif row1[0]=='alnCommand':
					self.alnCommand=row1[1]
				elif row1[0]=='alnTag':
					self.alnTag=int(row1[1])
			file1.close()
		else:
			self.alnFile,self.alnCommand,self.alnTag=self.minimap2()
		if self.alnTag!=0:
			self.alnFile,self.alnCommand,self.alnTag=self.minimap2()
		l='alnFile\t'+self.alnFile+"\nalnCommand\t"+self.alnCommand+"\nalnTag\t"+str(self.alnTag)+"\n"
		file1=open(self.selectRawReadsLog,'w')
		file1.writelines(l)
		file1.close()
		self.alnSortDepth=self.getDepth()

		self.readName,self.readFile,self.readLog=self.getReads()

	def getReads(self):
		readName=self.out+"/Genome.inputCtg.usefullReads.txt"
		readFile=self.out+"/Genome.inputCtg.usefullReads.fasta"
		readLog=self.out+"/Genome.inputCtg.usefullReads.log"
		if os.path.exists(readLog)==True and os.path.getsize(readLog)!=0:
			return readName,readFile,readLog
		else:
			file1=open(self.alnSortDepth,'r')
			depth=[]
			for row in file1:
				row1=row.rstrip().split('\t')
				depth.append(int(row1[-1]))
			file1.close()
			a=np.argmax(np.bincount(depth))
			
			a1=a*(2.0-self.filterDepth)
			a2=a*self.filterDepth
			readsname=[]
			file1=open(self.alnSortDepth,'r')
			loci={}
			for row in file1:
				row1=row.rstrip().split('\t')
				d=int(row1[-1])
				if d<=a2 or d>=a1:
					if loci.has_key(row1[0]):
						lo=loci[row1[0]]
					else:
						lo=[]
					lo.append(int(row1[1]))
					loci[row1[0]]=lo

			file1.close()
			
			file1=open(readName,'w')
			file2=open(readFile,'w')
			commendlines='samtools index -bc -@ '+self.thread+" "+self.alnFile
			os.system(commendlines)
			bamFile=AlignmentFile(self.alnFile,"rb")
		
			loci1={}
			for k,v in loci.items():
				v1=self.mergelist(v)
				loci[k]=v1
				print (k)
				print (len(v))
				print (len(v1))
				
				for i in v1:
					s,e=i
					for r in bamFile.fetch(k, s, e):
						l=r.query_name+"\n"
						file1.writelines(l)
			for k,v in self.inputSeqDict.items():
				print (k,v)
				for r in bamFile.fetch(k, 0, self.seedLen):
					l=r.query_name+"\n"
					file1.writelines(l)
				for r in bamFile.fetch(k, v-self.seedLen,v):
					l=r.query_name+"\n"
					file1.writelines(l)
					

			readlist=[]
			Chrlist=[]
			for r in bamFile:
				if r.reference_name not in Chrlist:
					Chrlist.append(r.reference_name)
					print (r.reference_name)
				gseq=self.readsDict[r.query_name]
				queryLen=len(gseq.seq)
				if r.is_unmapped==True:
					l=r.query_name+"\n"
					file1.writelines(l)
				else:
					if float(r.query_alignment_length)/queryLen<0.99:
						l=r.query_name+"\n"
						file1.writelines(l)
			file1.close()
			commondline='cat '+readName+" | sort | uniq >"+self.out+'/temp'
			os.system(commondline)
			commondline='mv '+self.out+'/temp '+readName
			os.system(commondline)
			file1=open(readName)
			for row in file1:
				row1=row.rstrip()
				gseq=self.readsDict[row1]
				l='>'+gseq.description+"\n"+gseq.seq+"\n"
				file2.writelines(l)
				
			file2.close()
			bamFile.close()
			file1=open(readLog,'w')
			l='depthThreshold(>)\t'+str(a)+'\nreadName\t'+readName+"\nreadFile\t"+readFile+"\n"
			file1.writelines(l)
			file1.close()
			return readName,readFile,readLog
		
	def mergelist(self,v):
		v1=[]
		it=[v[0],v[0]]
		for i in v[1:]:
			if i==it[1]+1:
				it[1]=i
			else:
				v1.append(it)
				it=[i,i]
		return v1

	def getDepth(self):
		alnSortName=self.out+"/Genome.inputCtg.mappedReads.sort.bam"
		commandline="samtools sort -@ "+self.thread+" -o "+alnSortName+" "+self.alnFile
		if os.path.exists(alnSortName)!=True or os.path.getsize(alnSortName)==0:
			output=os.system(commandline)
		alnSortDepth=self.out+"/Genome.inputCtg.mappedReads.sort.depth.txt"
		if os.path.exists(alnSortDepth)!=True or os.path.getsize(alnSortDepth)==0:
			commandline="samtools depth -aa  "+alnSortName+" >"+alnSortDepth
			output=os.system(commandline)
		#return alnSortName,alnSortDepth
		return alnSortDepth
			

	def minimap2(self):
		alnname=self.out+"/Genome.inputCtg.mappedReads.sort.bam"
		commandline="minimap2 -t "+self.thread+" -Y -ax asm20 "+self.inputSeq+" "+self.reads+" | samtools sort -@ "+self.thread+" -o "+alnname
		if os.path.exists(alnname)==True and os.path.getsize(alnname)!=0 and self.alnTag==0:
			return alnname,commandline,str(0)
		else:
			output=os.system(commandline)
			minimaptag=1
			if output!=0:
				while output!=0:
					output=os.system(commandline)
					minimaptag+=1
					if minimaptag>=3:
						print ("minimap2 cannot do proper alignment!!!")
						sys.exit()
			return alnname,commandline,str(output)	
