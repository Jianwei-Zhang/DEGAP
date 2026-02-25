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

class FindExtensionReads(object):
	def __init__(self,roundInput,lastRoundUsedReads,usedReads):
		self.roundInput=roundInput
		self.lastRoundUsedReads=lastRoundUsedReads
		self.usedReads=usedReads
		self.note=''
		self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
		self.log=self.roundInput.elongation.roundDir+"/extensionReads.log"
		if os.path.exists(self.extensionReads)==True and os.path.getsize(self.extensionReads)!=0:
			self.readlog()
		else:
			logfilet=open(self.log,'w')
			self.potentialExtensionReadsAln,self.minimap2Command,self.minimap2Output=self.minimap2()
			logLine='potentialExtensionReadsAln\t'+self.potentialExtensionReadsAln+"\nminimap2Command\t"+self.minimap2Command+"\nminimap2Output\t"+self.minimap2Output+"\n"
			logfilet.writelines(logLine)

			self.minimumExtensionReads()
			logLine='minimumThresholdReadsAln\t'+self.minimumThresholdReadsAln+"\n"
			logLine+='minimumThresholdReadsID\t'+';'.join(self.minimumThresholdReadsID)+"\n"
			logLine+='minimumThresholdExtensionReadsAln\t'+self.minimumThresholdExtensionReadsAln+"\n"
			logLine+='minimumThresholdExtensionReads\t'+self.minimumThresholdExtensionReads+"\n"
			logLine+='minimumThresholdExtensionReadsID\t'+';'.join(self.minimumThresholdExtensionReadsID)+"\n"
			logfilet.writelines(logLine)
			
			
			if self.note=='':
				self.selectMappingQuality=20
				self.selectAlignmentLength=3000
				self.selectNMAlignmentLengthratio=0.1
				
				self.selectReadsNum=0
				self.extensionReadsNum=0
				self.readsExtensionLength=1000
				self.extensionReadsEdge=10
			
				while self.extensionReadsNum<=0 and self.note=='':
					self.selectPotentialExtensionReadsAln=self.roundInput.elongation.roundDir+"/selectPotentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
					#self.selectPotentialExtensionReadsID=self.samFilter(self.potentialExtensionReadsAln,self.selectPotentialExtensionReadsAln)
					self.selectPotentialExtensionReadsID=self.samFilter(self.minimumThresholdExtensionReadsAln,self.selectPotentialExtensionReadsAln)
					self.selectReadsNum=len(self.selectPotentialExtensionReadsID)

					self.extensionReadsAln=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".bam"
					self.extensionReads=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
					self.extensionReadsID=self.extensionFinder(self.selectPotentialExtensionReadsAln,self.extensionReadsAln,self.extensionReads)
					self.extensionReadsNum=len(self.extensionReadsID)

					if self.extensionReadsNum==0:
						if self.selectMappingQuality==0:
							if self.readsExtensionLength<=10:
								if self.selectAlignmentLength<=500:
									self.extensionReadsEdge=self.extensionReadsEdge+10
									self.selectMappingQuality=20
									self.selectAlignmentLength=3000
									self.readsExtensionLength=1000
								else:
									self.selectAlignmentLength=self.selectAlignmentLength-100
							else:
								self.readsExtensionLength=self.readsExtensionLength-100
						else:
							self.selectMappingQuality=self.selectMappingQuality-10

				logLine='selectPotentialExtensionReadsAln\t'+self.selectPotentialExtensionReadsAln+"\n"
				logLine+='selectPotentialExtensionReadsID\t'+';'.join(self.selectPotentialExtensionReadsID)+"\n"
				logLine+='selectReadsNum\t'+str(self.selectReadsNum)+"\n"
				logLine+='extensionReadsAln\t'+self.extensionReadsAln+"\n"
				logLine+='extensionReads\t'+self.extensionReads+"\n"
				logLine+='extensionReadsID\t'+';'.join(self.extensionReadsID)+"\n"
				logLine+='extensionReadsNum\t'+str(self.extensionReadsNum)+"\n"
				logfilet.writelines(logLine)

			else:
				logLine='selectPotentialExtensionReadsAln\tNone\n'
				logLine+='selectPotentialExtensionReadsID\tNone\n'
				logLine+='selectReadsNum\t0\n'
				logLine+='extensionReadsAln\tNone\n'
				logLine+='extensionReads\tNone\n'
				logLine+='extensionReadsID\tNone\n'
				logLine+='extensionReadsNum\t0\n'
				logfilet.writelines(logLine)
			logLine='selectMappingQuality\t'+str(self.selectMappingQuality)+"\n"
			logLine+='selectAlignmentLength\t'+str(self.selectAlignmentLength)+"\n"
			logLine+='selectNMAlignmentLengthratio\t'+str(self.selectNMAlignmentLengthratio)+"\n"
			logLine+='readsExtensionLength\t'+str(self.readsExtensionLength)+"\n"
			logLine+='extensionReadsEdge\t'+str(self.extensionReadsEdge)+"\n"
			logLine+='note\t'+str(self.note)+"\n"
			logfilet.writelines(logLine)	

			logfilet.close()
	
	def minimumExtensionReads(self):
		self.selectMappingQuality=0
		self.selectAlignmentLength=500
		self.selectNMAlignmentLengthratio=0.1
		self.extensionReadsEdge=self.roundInput.elongation.base.edge
		self.readsExtensionLength=10

		self.minimumThresholdReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdReadsID=self.samFilter(self.potentialExtensionReadsAln,self.minimumThresholdReadsAln)


		self.minimumThresholdExtensionReadsAln=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		self.minimumThresholdExtensionReads=self.roundInput.elongation.roundDir+"/minimumThresholdExtensionReads."+self.roundInput.elongation.base.tag+".fa"
		self.minimumThresholdExtensionReadsID=self.extensionFinder(self.minimumThresholdReadsAln,self.minimumThresholdExtensionReadsAln,self.minimumThresholdExtensionReads)

		if len(self.minimumThresholdExtensionReadsID)==0:
			self.note='noExtensionReadsFoundAtMinimumThreshold'
		else:
			self.note=''

	def extensionFinder(self,inputAln,outputAln,outputSeq):
		readslist=[]
		outputSeqFile=open(outputSeq,'w')
		inputAlnFile=AlignmentFile(inputAln,"rb",check_sq=False)
		outputAlnFile=AlignmentFile(outputAln,"wb", template=inputAlnFile)
		#print inputAln,outputAln
		readsExtensionLength=self.readsExtensionLength
		extensionReadsEdge=self.extensionReadsEdge

		for r in inputAlnFile:
			queryID=r.query_name
			query=self.roundInput.elongation.base.readsDict[queryID]
			queryDistance,refDistance,extensionLength,extensionReadsSeq=self.calculateBoundDistance(queryID,query,r)
			if self.readsExtensionLength<=extensionLength:
				if queryDistance<=self.extensionReadsEdge and refDistance<=self.extensionReadsEdge:
					if queryID not in self.usedReads :#or queryID in self.lastRoundUsedReads:
						if queryID not in readslist:
							l='>'+queryID+"\n"+query.seq+"\n"
							outputSeqFile.writelines(l)
							readslist.append(queryID)
					outputAlnFile.write(r)
		
		outputSeqFile.close()
		inputAlnFile.close()
		outputAlnFile.close()
		return readslist
	
	def calculateBoundDistance(self,queryID,query,r):
		if r.is_reverse:
			queryseq=query.seq.reverse_complement()
		else:
			queryseq=query.seq
		
		qs,qe,rs,re=self.findAlnPosition(queryseq,self.roundInput.inputSeedSequence,r)
		if self.roundInput.elongation.base.flag=='left':
			queryDistance=qs
			extensionReadsSeq=queryseq
			refDistance=len(self.roundInput.inputSeedSequence.seq)-re
			extensionLength=len(queryseq)-qe
		else:
			queryDistance=len(queryseq)-qe
			extensionReadsSeq=queryseq
			refDistance=rs
			extensionLength=qs
		if qs<0 or len(queryseq)<qe or qe<0 or rs<0 or re<0 :#or len(self.roundInput.inputSeedSequence.seq)<re:
			print ('wrong seqend',qs,qe,rs,re)
			print (queryDistance,refDistance,extensionLength,extensionReadsSeq)
			print (qs<0)
			print (len(queryseq)<qe)
			print (qe<0)
			print (rs<0)
			print (re<0)
			print (len(self.roundInput.inputSeedSequence.seq)<re)
			print (len(self.roundInput.inputSeedSequence.seq),re)
			sys.exit()
		return queryDistance,refDistance,extensionLength,extensionReadsSeq
	
	def findAlnPosition(self,queryseq,refseq,r):
		ct=r.cigartuples
		scpstart=ct[0]
		if scpstart[0]==4:
			scps=scpstart[1]
		else:
			scps=0
	
		scpend=ct[-1]
		if scpend[0]==4:
			scpe=scpend[1]
		else:
			scpe=0

		if queryseq==r.query_sequence or r.query_alignment_sequence==None:
			qs=r.query_alignment_start
			qe=r.query_alignment_end
		else:
			qaln=r.query_alignment_sequence
			qs=str(queryseq).index(qaln)
			qe=qs+len(r.query_alignment_sequence)
		rs=r.reference_start
		re=r.reference_end
		if qs<0:
			print ('wrong qs',qs)
		if qe>len(queryseq):
			print ('wrong qe',qe,len(queryseq))
		return qs,qe,rs,re

		

	def samFilter(self,inputAln,outputAln):
		samFile=AlignmentFile(inputAln,"rb",check_sq=False)
		outputBamFile=AlignmentFile(outputAln,"wb", template=samFile)
		
		readslist=[]
		print (inputAln)
		sMQ=self.selectMappingQuality
		sAlignmentLength=self.selectAlignmentLength
		sNMAlignmentLengthr=self.selectNMAlignmentLengthratio

		for r in samFile:
			if r.is_unmapped==False and r.mapping_quality>=sMQ:
				for i in r.tags:
					if i[0]=='NM':
						NM=i[1]
				AlignmentLength=len(r.query_alignment_sequence)
				if  AlignmentLength>=sAlignmentLength and float(NM)/AlignmentLength<=sNMAlignmentLengthr:
					outputBamFile.write(r)
					if r.query_name not in readslist:
						readslist.append(r.query_name)
		samFile.close()
		outputBamFile.close()
		return readslist
		
	def minimap2(self):
		alnname=self.roundInput.elongation.roundDir+"/potentialExtensionReads."+self.roundInput.elongation.base.tag+".bam"
		alnname1=self.roundInput.elongation.roundDir+"/extensionReads."+self.roundInput.elongation.base.tag+".fa"
		commandline="minimap2 -t "+self.roundInput.elongation.base.thread+" -Y -ax asm20 "+self.roundInput.inputSeq+" "+self.roundInput.elongation.base.reads+" | samtools view -bS >"+alnname
		if os.path.exists(alnname)==True and os.path.getsize(alnname)!=0 and os.path.exists(alnname1)==True and os.path.getsize(alnname1)!=0:
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
			


	def readlog(self):
		logfilet=open(self.log,'r')
		for row in logfilet:
			row1=row.rstrip().split('\t')
			if row1[0]=='potentialExtensionReadsAln':
				self.potentialExtensionReadsAln=row1[1]
			elif row1[0]=='minimap2Command':
				self.minimap2Command=row1[1]
			elif row1[0]=='minimap2Output':
				self.minimap2Output=row1[1]
			elif row1[0]=='minimumThresholdReadsAln':
				self.minimumThresholdReadsAln=row1[1]
			elif row1[0]=='minimumThresholdReadsID':
				self.minimumThresholdReadsID=row1[1].split(';')
			elif row1[0]=='minimumThresholdExtensionReadsAln':
				self.minimumThresholdExtensionReadsAln=row1[1]
			elif row1[0]=='minimumThresholdExtensionReads':
				self.minimumThresholdExtensionReads=row1[1]
			elif row1[0]=='minimumThresholdExtensionReadsID':
				self.minimumThresholdExtensionReadsID=row1[1].split(';')
			elif row1[0]=='selectPotentialExtensionReadsAln':
				self.selectPotentialExtensionReadsAln=row1[1]

			elif row1[0]=='selectPotentialExtensionReadsID':
				if row1[1]!='None':
					self.selectPotentialExtensionReadsID=row1[1].split(';')
				else:
					self.selectPotentialExtensionReadsID=[]
			elif row1[0]=='selectReadsNum':
				self.selectReadsNum=int(row1[1])
			elif row1[0]=='extensionReadsAln':
				self.extensionReadsAln=row1[1]
			elif row1[0]=='extensionReads':
				self.extensionReads=row1[1]
			elif row1[0]=='extensionReadsID':
				if row1[1]!='None':
					self.extensionReadsID=row1[1].split(';')
				else:
					self.extensionReadsID=[]
			elif row1[0]=='extensionReadsNum':
				self.extensionReadsNum=int(row1[1])
			elif row1[0]=='selectMappingQuality':
				self.selectMappingQuality=int(row1[1])
			elif row1[0]=='selectAlignmentLength':
				self.selectAlignmentLength=int(row1[1])
			elif row1[0]=='selectNMAlignmentLengthratio':
				self.selectNMAlignmentLengthratio=float(row1[1])
			elif row1[0]=='readsExtensionLength':
				self.readsExtensionLength=int(row1[1])
			elif row1[0]=='extensionReadsEdge':
				self.extensionReadsEdge=int(row1[1])
			elif row1[0]=='note':
				if len(row1)==1:
					self.note=''
				else:
					self.note=row1[1]
		logfilet.close()
