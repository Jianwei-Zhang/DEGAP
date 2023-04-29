import re
import os
import sys
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO
import GapFillerClass
from GapFillerClass import GapFillerClass,mummer

class GapFiller(object):
	def __init__(self,parameterlist):
		self.mode,self.remove,self.thread,self.reads,self.out,self.seqLeft,self.seqRight,self.flag,self.edge,self.filterDepth,self.readsDict,self.maxReadsLen,self.seedLen=parameterlist
		out=self.out
		self.log=out+"/process.log"
		if not os.path.exists(self.log):
			open(self.log,'w').close()
		self.summary=out+"/process.summary"
		if not os.path.exists(self.summary):
			open(self.summary,'w').close()
		self.agp=out+"/process.agp"
		self.usedReads=[]
		if not os.path.exists(self.agp):
			open(self.agp,'w').close()
		self.name=out.split('/')[-1]
		self.outfile=out+"/process"
		if not os.path.exists(self.outfile):
			os.makedirs(self.outfile)

		if self.flag=='left':
			self.initialSeq=self.seqLeft
			self.terminalSeq=self.seqRight
			self.tag='left'
		else:
			self.initialSeq=self.seqRight
			self.terminalSeq=self.seqLeft
			self.tag='right'

		self.Elongation=Elongation(self)
	
class Elongation(object):
	def __init__(self,base):
		self.base=base
		self.roundNum=1
		self.usedReads=[]
		self.extensionLen=0
		self.endSignal=False

		logfile=open(self.base.log,'w')
		summaryfile=open(self.base.summary,'w')
		while self.endSignal==False:
			self.ElongationInit(logfile)
			self.lastRoundUsedReads=self.ElongateSeq(logfile,summaryfile,self.lastRoundUsedReads)
			if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
				self.extensionLen=self.extensionLen+self.roundResult.ExtensionContigs.extensionLength
			if self.roundNum==1:
				if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
					atgInitial=self.roundResult.ExtensionContigs.selectExtensionContigsAln[0][0]
				else:
					atgInitial=''
			self.removeFile()
			self.roundNum+=1
		summaryfile.close()

		self.roundNum=self.roundNum-1
		self.finalSeq=self.base.out+'/'+self.base.name+'.final.fa'
		if os.path.exists(self.roundDir+"/linkedSequence.fasta")==True and  os.path.getsize(self.roundDir+"/linkedSequence.fasta")!=0:
			fileofs=open(self.finalSeq,'w')
			l='>'+self.base.name+"\n"
			fileofs.writelines(l)
			for gseq in SeqIO.parse(self.roundResult.roundOutput.linkedSequence,'fasta'):
				l0=gseq.description
				l1=l0.split('\t')
				l=gseq.seq+'\n'
				fileofs.writelines(l)
				for l11 in l1:
					l2=l11.split(':')
					if l2[0]=="Aln":
						atgTerminal=l2[1]
			fileofs.close()
		elif self.roundResult.ExtensionReads.note!='' or 'No extension contigs or reads found' in self.roundResult.ExtensionContigs.selectContigNote:
			fileofs=open(self.finalSeq,'w')
			l='>'+self.base.name+"_noExtensionContigsorReads\n"
			fileofs.writelines(l)
			atgTerminal=''
			for gseq in SeqIO.parse(self.roundInputSeq,'fasta'):
				l=gseq.seq+'\n'
				fileofs.writelines(l)
			fileofs.close()
		else:
			fileofs=open(self.finalSeq,'w')
			l='>'+self.base.name+"_noNewExtensionReads\n"
			fileofs.writelines(l)
			atgTerminal=''
			for gseq in SeqIO.parse(self.roundResult.roundOutput.outputSequence,'fasta'):
				l=gseq.seq+'\n'
				fileofs.writelines(l)
			fileofs.close()
	
		for fgseq in SeqIO.parse(self.finalSeq,'fasta'):
			l='Final ExtensionSequence: '+fgseq.id+"\n"
			l=l+'Final EXtendFile: '+self.finalSeq+"\n"
			logfile.writelines(l)
		logfile.close()
		
		agpfile=open(self.base.agp,'w')
		for gseq in SeqIO.parse(self.base.initialSeq,'fasta'):
			initialSequence=gseq
		if atgInitial!='':
			if self.base.flag=='left':
				atgInitial1=atgInitial.split('\t')
				st=len(initialSequence.seq)-(int(atgInitial1[7])-int(atgInitial1[1]))
				if atgTerminal=='':
					l=self.base.name+"\t1"+"\t"+str(st)+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq.seq))+"\t2\tw\t"+fgseq.id+"\t"+str(st+1)+"\t"+str(len(fgseq.seq))+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=None
				else:
					l=self.base.name+"\t1"+"\t"+str(st)+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					atgTerminal1=atgTerminal.split(';')
					eT=len(fgseq.seq)-(int(atgTerminal1[7])-int(atgTerminal1[4]))
					l=self.base.name+"\t"+str(st+1)+"\t"+str(eT)+"\t2\tw\t"+fgseq.id+"\t"+str(st+1)+"\t"+str(eT)+"\t+\n"
					agpfile.writelines(l)
					Terminalname1=atgTerminal1[-2]
					l=self.base.name+"\t"+str(int(eT)+1)+"\t"+str(len(fgseq.seq))+"\t3\tw\t"+Terminalname1+"\t"+str(int(atgTerminal1[1])+1)+"\t"+atgTerminal1[7]+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=Terminalname1
	
			else:
				atgInitial1=atgInitial.split('\t')
				st=len(fgseq.seq)-len(initialSequence.seq)+int(atgInitial1[0])-1
				if atgTerminal=='':
					atgInitial1=atgInitial.split('\t')
					l=self.base.name+"\t1\t"+str(st)+"\t1\tw\t"+fgseq.id+"\t1\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq))+"\t2\tw\t"+initialSequence.id+"\t"+atgInitial1[0]+"\t"+str(len(initialSequence.seq))+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=None
					
				else:
					atgTerminal1=atgTerminal.split(';')
					atgInitial1=atgInitial.split('\t')
					sT=int(atgTerminal1[1])
					l=self.base.name+"\t1\t"+str(sT)+"\t1\tw\t"+atgTerminal1[-2]+"\t1\t"+atgTerminal1[1]+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(sT+1)+"\t"+str(st)+"\t2\tw\t"+fgseq.id+"\t"+str(sT+1)+"\t"+str(st)+"\t+\n"
					agpfile.writelines(l)
					l=self.base.name+"\t"+str(st+1)+"\t"+str(len(fgseq))+"\t3\tw\t"+initialSequence.id+"\t1\t"+atgTerminal1[1]+"\t+\n"
					agpfile.writelines(l)
					self.TerminalSeq=atgTerminal1[-2]
		else:
			l=self.base.name+"\t1"+"\t"+str(len(initialSequence.seq))+"\t1\tw\t"+initialSequence.id+"\t1\t"+str(len(initialSequence.seq))+"\t+\n"
			agpfile.writelines(l)
		agpfile.close()
		self.removeFile()
				
	def removeFile(self):
		if os.path.exists(self.roundDir+"/hifiasm")==True or os.path.exists(self.roundDir+"/potentialExtensionReads.left.sam")==True:
			if self.base.remove==2 or self.base.remove==1:
				commondline='rm '+self.roundDir+"/*.bam"
				os.system(commondline)
				commondline='rm -rf '+self.roundDir+"/hifiasm"
				os.system(commondline)
		if self.endSignal==True:
			if os.path.exists(self.base.outfile)==True:
				if self.base.remove==1:
					commondline='rm -rf '+self.base.outfile
					os.system(commondline)
		
	

	def ElongateSeq(self,logfile,summaryfile,lastRoundUsedReads):
		self.lastRoundUsedReads=lastRoundUsedReads
		roundLog=open(self.roundLog,'w')
		roundSummary=open(self.roundSummary,'w')
		self.roundResult=GapFillerClass(self)
		logLine,summeryLine=self.writelog()
		print (logLine)
		roundLog.writelines(logLine)
		roundSummary.writelines(summeryLine)
		logfile.writelines(logLine)
		summaryfile.writelines(summeryLine)
		roundLog.close()
		roundSummary.close()
		if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
			return self.roundResult.roundOutput.ExtensionUsedReads
		else:
			return []
	
	def writelog(self):
		logLine='\n\n*****************\n\n'
		logLine+='\toutputPath: '+str(self.roundDir)+"\n"
		logLine+="\tseedSequenceLength: "+str(self.base.seedLen)+"\n"
		logLine+="\tinitialSequenceFile: "+str(self.base.initialSeq)+"\n"
		for gseq in SeqIO.parse(self.roundInputSeq,'fasta'):
			inputSeq=gseq
		logLine+="\t\tinitialSeqnenceID: "+str(inputSeq.id)+"\n\t\tinitialSeqnenceLength: "+str(len(inputSeq.seq))+"\n"
		logLine+="\tterminalSequenceFile: "+str(self.base.terminalSeq)+"\n"
		logLine+="\tseedSequenceFile: "+str(self.roundResult.roundInput.inputSeq)+"\n"
		logLine+="\t\tseedSeqnenceID: "+str(self.roundResult.roundInput.inputSeedSequence.id)+"\n\t\tseedSeqnenceLength: "+str(len(self.roundResult.roundInput.inputSeedSequence.seq))+"\n\n"
		if self.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in self.roundResult.ExtensionContigs.selectContigNote:
			logLine+='minimap2Commond: '+str(self.roundResult.ExtensionReads.minimap2Command)+"\n"
			logLine+='\textensionReads: \n\tselectReadsNum: '+str(self.roundResult.ExtensionReads.selectReadsNum)+"\n"
			logLine+="\t\tselectReadsAln: "+str(self.roundResult.ExtensionReads.selectPotentialExtensionReadsAln)+"\n"
			logLine+="\t\tselectMappingQuality: "+str(self.roundResult.ExtensionReads.selectMappingQuality)+"\n"
			logLine+="\t\tselectAlignmentLength: "+str(self.roundResult.ExtensionReads.selectAlignmentLength)+"\n"
			logLine+="\t\tselectNMAlignmentLengthratio: "+str(self.roundResult.ExtensionReads.selectNMAlignmentLengthratio)+"\n\n"
			logLine+="\textensionReadsNum: "+str(self.roundResult.ExtensionReads.extensionReadsNum)+"\n"
			logLine+="\t\textensionReadsAln: "+str(self.roundResult.ExtensionReads.extensionReadsAln)+"\n"
			logLine+="\t\textensionReadsFile: "+str(self.roundResult.ExtensionReads.extensionReads)+"\n"
			logLine+="\t\textensionReadsMinimumExtensionLength: "+str(self.roundResult.ExtensionReads.readsExtensionLength)+"\n"
			logLine+="\t\textensionReadsMaximumEdge: "+str(self.roundResult.ExtensionReads.extensionReadsEdge)+"\n\n"

			logLine+="\textensionSequnece: "+str(self.roundResult.ExtensionContigs.extensionContigs)+"\n"
			logLine+="\textensionSequneceNote: "+str(self.roundResult.ExtensionContigs.extensionSeqNote)+"\n"+str(self.roundResult.ExtensionContigs.selectContigNote)+"\n"
			logLine+="\t\textensionSequneceIdentity: "+str(self.roundResult.ExtensionContigs.selectContigIdentity)+"\n"
			logLine+="\t\textensionSequneceMinimumExtensionLength: "+str(self.roundResult.ExtensionContigs.selectContigAlnLength)+"\n"
			logLine+="\t\textensionSequneceMaximumEdge: "+str(self.roundResult.ExtensionContigs.selectContigDistance)+"\n"
			logLine+="\t\textensionSequneceAlnMerge: "+str(self.roundResult.ExtensionContigs.contigAlnMerge)+"\n"
			logLine+="\t\textensionSequneceAlnMergeIdentity: "+str(self.roundResult.ExtensionContigs.contigAlnMergeIdentity)+"\n"
			logLine+="\textensionedSeedSequenceFile: "+str(self.roundResult.ExtensionContigs.extensionSequence)+"\n"
			logLine+="\textensionedLength: "+str(self.roundResult.ExtensionContigs.extensionLength)+"\n"
			logLine+="\t\textensionContigOrReadsID:\n\t\t\t"+'\n\t\t\t'.join(self.roundResult.ExtensionContigs.extensionContigID)+"\n"
			newReads,note=self.updateUsedReads()
			logLine+=note
			logLine+="\t\tusedReadsNum: "+str(len(self.roundResult.roundOutput.ExtensionUsedReads))+"\n"
			logLine+="\t\tusedReads:\n\t\t\t"+"\n\t\t\t".join(self.roundResult.roundOutput.ExtensionUsedReads)+"\n"
			logLine+="\t\tusedNewReads:\n\t\t\t"+"\n\t\t\t".join(newReads)+"\n\n"
			logLine+="\toutputFile: "+str(self.roundResult.roundOutput.outputSequence)+"\n"
			logLine+="\t\toutputSequenceLength: "+str(self.roundResult.roundOutput.totalOutputSequenceLength)+"\n\n"
			
			if self.roundResult.roundOutput.linkedSequenceNote!='':
				logLine+='\tGAP can be closed!\n'+str(self.roundResult.roundOutput.linkedSequenceNote)+"\nLinkedSequence File: "+str(self.roundResult.roundOutput.linkedSequence)+"\n"
				logLine+='Endloop!\t'+str(self.roundResult.roundOutput.linkedSequence)+"\n"
				for gseq in SeqIO.parse(self.roundResult.roundOutput.linkedSequence,'fasta'):
					l0=gseq.description
					l1=l0.split('\t')
					for l11 in l1:
						l2=l11.split(':')
						if l2[0]=="Aln":
							closectg=l2[1].split(';')
							logLine+="Linked ctg:\t"+closectg[-2]+"\n"
							logLine+='\t'.join(closectg)+"\n"
						self.endSignal=True
			else:
				logLine+='\tGAP still not closed!\n'
			logLine+='\n\n*****************\n\n'
			summeryLine='round'+str(self.roundNum)+"\t"+str(len(inputSeq.seq))+"\t"+str(self.roundResult.roundOutput.totalOutputSequenceLength)+"\t"+str(self.roundResult.ExtensionContigs.extensionLength)+"\t"+"-ovl-".join(self.roundResult.ExtensionContigs.extensionContigID)+"\t"+str(len(newReads))+"\t"+";".join(newReads)+"\t"+str(len(self.roundResult.roundOutput.ExtensionUsedReads))+"\t"+";".join(self.roundResult.roundOutput.ExtensionUsedReads)+"\n"
			return logLine,summeryLine
		else:
			logLine+='No ExtensionReads or  ExtensionContig Found\n'
			logLine+="Endloop!\t"+self.roundResult.ExtensionReads.note+"\n"
			self.endSignal=True
			summeryLine=''
			return logLine,summeryLine

	def updateUsedReads(self):
		newReads=[]
		sameWithLastRound=[]
		note=''
		for i1 in self.roundResult.roundOutput.ExtensionUsedReads:
			if i1 not in self.usedReads:
				self.usedReads.append(i1)
				newReads.append(i1)
			if i1 in self.lastRoundUsedReads:
				sameWithLastRound.append(i1)
		if len(newReads)==0 :
			note='No New ExtensionReads Found,\t'
			if len(sameWithLastRound)==0:
				note+='Not same ExtensionReads with last round ExtensionReads, end up a loop!!\n'
				self.endSignal=True
			else:
				note+='Same ExtensionReads with the last round,continune the loop!\n'
				if self.roundResult.ExtensionContigs.extensionLength==0:
					note+='However,ExtensionLength==0,end up a loop!!\n'
					self.endSignal=True
		return newReads,note

	
	def ElongationInit(self,logfile):
		l='\n\n*****************\n\nExtensionRound '+str(self.roundNum)+'\n'
		logfile.writelines(l)
		print (l)
		self.roundDir=self.base.outfile+'/round'+str(self.roundNum)
		if not os.path.exists(self.roundDir):
			os.makedirs(self.roundDir)
		self.lastRoundDir=self.base.outfile+'/round'+str(self.roundNum-1)
		self.roundLog=self.base.outfile+'/round'+str(self.roundNum)+"/log"
		self.roundSummary=self.base.outfile+'/round'+str(self.roundNum)+"/summary"

		if self.roundNum!=1:
			self.roundInputSeq=self.lastRoundDir+'/outputExtensionSequence.fasta'
		else:
			self.roundInputSeq=self.base.initialSeq
			self.lastRoundUsedReads=[]

