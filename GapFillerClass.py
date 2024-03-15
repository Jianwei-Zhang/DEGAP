import re
import Bio
import os
import sys
import re
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO
import FindExtensionReads
from FindExtensionReads import FindExtensionReads
import FindExtensionContigs
from FindExtensionContigs import FindExtensionContigs

class InputSequence(object):
	def __init__(self,Elongation):
		self.elongation=Elongation
		self.inputSeq=self.elongation.roundDir+"/inputCutSequence.fasta"
		filet=open(self.inputSeq,'w')
		for gseq in SeqIO.parse(self.elongation.roundInputSeq,'fasta'):
			length=int(self.elongation.base.seedLen)
			if len(gseq.seq)<length:
				l='>'+gseq.id+"_all\n"+gseq.seq+"\n"
			else:
				if self.elongation.base.flag=='left':
					l='>'+gseq.id+"_"+str(length)+"\n"+gseq.seq[-length:]+"\n"
				else:
					l='>'+gseq.id+"_"+str(length)+"\n"+gseq.seq[:length]+"\n"
			filet.writelines(l)
		filet.close()
		for gseq in SeqIO.parse(self.inputSeq,'fasta'):
			self.inputSeedSequence=gseq

def mummer(seq1,seq2,name):
	mr=name
	commendline="nucmer -c 90 -l 40 -p "+mr+" "+seq1+" "+seq2	
	os.system(commendline)
	commendline="delta-filter -m "+mr+".delta > "+mr+".delta.filter"
	os.system(commendline)
	commendline="show-coords -TrHcl "+mr+".delta.filter > "+mr+".delta.filter.coords"
	os.system(commendline)
	mr=mr+".delta.filter.coords"
	return mr

class OutputSequence(object):
	def __init__(self,ExtensionContigs,Elongation):
		self.extensionContigs=ExtensionContigs
		self.Elongation=Elongation
		self.outputLog=self.Elongation.roundDir+"/outputExtensionSequence.log"
		self.outputSequence=self.Elongation.roundDir+"/outputExtensionSequence.fasta"
		if os.path.exists(self.outputSequence)==True and os.path.getsize(self.outputSequence)!=0:
                        self.readlog()
		else:
			logfile=open(self.outputLog,'w')
			self.totalOutputSequenceLength=self.extensionSequence()
			logLine="outputSequence\t"+self.outputSequence+"\n"
			logLine+="totalOutputSequenceLength\t"+str(self.totalOutputSequenceLength)+"\n"
			self.ExtensionUsedReads=self.checkExtensionReads()
			logLine+="ExtensionUsedReads\t"+";".join(self.ExtensionUsedReads)+"\n"
			self.linkedSequence,self.linkedSequenceNote,self.linkedSequenceAln=self.linkgap()
			logLine+="linkedSequence\t"+self.linkedSequence+"\n"
			for i in self.linkedSequenceNote.rstrip().split('\n'):
				logLine+="linkedSequenceNote\t"+i+"\n"
			logLine+="linkedSequenceAln\t"+self.linkedSequenceAln+"\n"
			logfile.writelines(logLine)
			logfile.close()
	def readlog(self):
		logfile=open(self.outputLog,'r')
		self.linkedSequenceNote=''
		for row in logfile:
			row1=row.rstrip().split('\t')
			if row1[0]=='outputSequence':
				self.outputSequence=row1[1]
			elif row1[0]=="totalOutputSequenceLength":
				self.totalOutputSequenceLength=int(row1[1])
			elif row1[0]=="ExtensionUsedReads":
				self.ExtensionUsedReads=row1[1].split(';')
			elif row1[0]=="linkedSequenceNote":
				self.linkedSequenceNote=self.linkedSequenceNote+"\t".join(row1[1:])
			elif row1[0]=='linkedSequence':
				self.linkedSequence=row1[1]
			elif row1[0]=="linkedSequenceAln":
				self.linkedSequenceAln="\t".join(row1[1:])
		logfile.close()
	
	def linkgap(self):
		mummern=self.Elongation.roundDir+"/extensionSequence.linkgap.aln."+self.Elongation.base.tag+".mummer"
		mummerout=mummern+".delta.filter.coords"
		if os.path.exists(mummerout)!=True:
			mummerout=mummer(self.Elongation.base.terminalSeq,self.outputSequence,mummern)
		n1=self.Elongation.roundDir+"/linkedSequence.fasta"
		note=''
		f1=open(mummerout,'r')
		lout1=[]
		for row in f1:
			row1=row.rstrip().split('\t')
			len1=int(row1[4])
			if len1>=10000:
				lout1.append(row)

		f1.close()

		lout2=[]
		flag=self.Elongation.base.flag
		edg=self.Elongation.base.edge
		for row in lout1:
			row1=row.rstrip().split('\t')
			a,b,c,d=int(row1[0]),int(row1[1]),int(row1[2]),int(row1[3])
			ctgmame=row1[-1]
			if flag=='left':
				if d>=int(row1[8])-edg and a<=edg:
					lout2.append(row)
			else:
				if b>=int(row1[7])-edg and c<=edg:
					lout2.append(row)
				
		lout=''
		idty=0
		for row in lout2:
			row1=row.rstrip().split('\t')
			e=float(row1[6])
			if lout==[]:
				lout=row
				idty=e
			else:
				if e>idty:
					idty=e
					lout=row
		
		if lout!='':
			note=note+'ExtensionSequence can close the GAP!\n'
			self.Elongation.endSignal=True
			row1=lout.rstrip().split('\t')
			a,b,c,d=int(row1[0]),int(row1[1]),int(row1[2]),int(row1[3])
			for gseq in SeqIO.parse(self.outputSequence,'fasta'):
				Seq2=gseq.seq
			for gseq in SeqIO.parse(self.Elongation.base.terminalSeq,'fasta'):
				if gseq.id==row1[-2]:
					Seq1=gseq.seq
			for gseq in SeqIO.parse(self.Elongation.base.initialSeq,'fasta'):
				Seq0=gseq.seq
			if flag!='left':
				seqfin=Seq1[:b]+Seq2[d:]
				note=note+'\tGAP Length: '+str(self.Elongation.extensionLen+self.extensionContigs.extensionLength-int(row1[5]))+"\n\tLinked Sequence Length: "+str(len(seqfin))+"\n"
				l='>ExtensionSequence'+"\tgap:"+str(self.Elongation.extensionLen+self.extensionContigs.extensionLength-int(row1[5]))+"\tlen:"+str(len(seqfin))+"\tRound:round"+self.Elongation.roundDir.split('round')[-1]+"\tAln:"+';'.join(row1)+"\n"+seqfin+"\n"
			else:
				seqfin=Seq2+Seq1[b:]
				note=note+'\tGAP Length: '+str(self.Elongation.extensionLen+self.extensionContigs.extensionLength-int(row1[5]))+"\n\tLinked Sequence Length: "+str(len(seqfin))+"\n"
				l='>ExtensionSequence'+"\tgap:"+str(self.Elongation.extensionLen+self.extensionContigs.extensionLength-int(row1[5]))+"\tlen:"+str(len(seqfin))+"\tRound:round"+self.Elongation.roundDir.split('round')[-1]+"\tAln:"+';'.join(row1)+"\n"+seqfin+"\n"
			ft=open(n1,'w')
			ft.writelines(l)
			ft.close()
		return n1,note,lout

	def checkExtensionReads(self):
		IDlist=self.extensionContigs.extensionContigID
		readslist=[]
		for i in IDlist:
			if i in self.extensionContigs.hifiasmResultDict:
				readslist=readslist+self.extensionContigs.hifiasmResultDict[i]
			else:
				readslist.append(i)
		return readslist

	def extensionSequence(self):
		ft=open(self.outputSequence,'w')
		seedlen=int(self.Elongation.base.seedLen)
		i=0
		if 'No extension contigs or reads found' not in self.extensionContigs.selectContigNote:
			for gseq in SeqIO.parse(self.extensionContigs.extensionSequence,'fasta'):
				for gseq1 in SeqIO.parse(self.Elongation.roundInputSeq,'fasta'):
					if self.Elongation.base.flag=='left':
						Seq1=gseq1.seq[:-seedlen]
						seq2=gseq1.seq[:-seedlen]+gseq1.seq[-seedlen:]
						if not gseq1.seq==seq2:
							print (gseq1.seq==seq2)
							print ('wrong link')
							sys.exit()
						Seq2=Seq1+gseq.seq
						l='>'+gseq.id+'\n'+Seq2+"\n"
						ft.writelines(l)
					else:
						Seq1=gseq1.seq[seedlen:]
						seq2=gseq1.seq[:seedlen]+gseq1.seq[seedlen:]
						if not gseq1.seq==seq2:
							print (gseq1.seq==seq2)
							print ('wrong link')
							sys.exit()
						Seq2=gseq.seq+Seq1
						l='>'+gseq.id+'\n'+Seq2+"\n"
						ft.writelines(l)
						i+=1
			if i>=2:
				for gseq in SeqIO.parse(self.extensionContigs.extensionSequence,'fasta'):
					for gseq1 in SeqIO.parse(self.Elongation.roundInputSeq,'fasta'):
						print (gseq.id,gseq1.id)
				sys.exit()
			ft.close()
		else:
			for gseq1 in SeqIO.parse(self.Elongation.roundInputSeq,'fasta'):
				l='>'+gseq.id+'\n'+gseq.seq+"\n"
				Seq2=gseq.seq
				ft.writelines(l)
		ft.close()
		return len(Seq2)



class GapFillerClass(object):
	def __init__(self,Elongation):
		self.usedReads=Elongation.usedReads
		self.lastRoundUsedReads=Elongation.lastRoundUsedReads
		self.Elongation=Elongation
		self.roundInput=InputSequence(self.Elongation)
		self.ExtensionReads=FindExtensionReads(self.roundInput,self.lastRoundUsedReads,self.usedReads)
		if self.ExtensionReads.note=='':
			self.ExtensionContigs=FindExtensionContigs(self.ExtensionReads)
			if 'No extension contigs or reads found' not in self.ExtensionContigs.selectContigNote and "Reach the maximum Length" not in self.ExtensionContigs.selectContigNote:
				self.roundOutput=OutputSequence(self.ExtensionContigs,self.Elongation)
