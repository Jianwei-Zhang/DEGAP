import re
import os
import sys
import getopt
import pysam
from pysam import AlignmentFile
import Bio
from Bio import SeqIO
import GapFiller
from GapFiller import GapFiller
import N50Check
from N50Check import N50Check

class CtgLinker(object):
	def __init__(self,parameterlist):
		self.mode,self.remove,self.thread,self.reads,self.out,self.genomeSeq,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict,self.maxReadsLen,self.seedLen=parameterlist
		self.ctgSeq=self.out+"/Genome.inputCtg.fa"
		self.ctgSeqLog=self.out+"/Genome.inputCtg.log"
		if os.path.exists(self.ctgSeqLog)!=True or os.path.getsize(self.ctgSeqLog)==0:
			N50Check([self.genomeSeq,self.seedLen,self.ctgSeq,self.ctgSeqLog])
		self.projectSeq=self.out+"/DG.project.fa"
		file1= open(self.projectSeq, 'w').close()

		self.project=self.out+"/project"
		if not os.path.exists(self.project):
			os.makedirs(self.project)

		self.unplaced=self.out+"/unplaced.txt"
		file1=open(self.unplaced,'w')
		i=0
		for gseq in SeqIO.parse(self.ctgSeq,'fasta'):
			l=gseq.description+"\n"
			file1.writelines(l)	
			i+=1
		self.unplacednum=i
		file1.close()

		self.used=self.out+"/used.txt"
		file1=open(self.used,'w')
		l='# "*" means sequences still can elongate\n# "**" means sequences can no longer elongate\n#DG\tright-end\tleft-end\ttags\n'
		file1.writelines(l)
		file1.close()

		self.agppwd=self.out+"/DG.agp.path.txt"
		file1=open(self.agppwd,'w').close()

		self.agp=self.out+"/DG.scaffold.agp"
		file1=open(self.agp,'w').close()

		self.scaffoldSeq=self.out+"/DG.scaffold.fa"
		file1=open(self.scaffoldSeq,'w').close()

		self.CtgLinkerPipline()
		self.scaffoldSeqLog=self.out+"/DG.scaffold.log"
		N50Check([self.scaffoldSeq,self.scaffoldSeqLog])

	def CtgLinkerPipline(self):
		self.scaffoldID=1
		self.projectID=1
		self.Lastround=False
		while self.unplacednum >0 or self.Lastround==False:
			self.roundInit()
			if os.path.exists(self.projectfasta)==True and os.path.getsize(self.projectfasta)!=0:
				self.skipRound()
			else:
				self.elongationDG()
			self.projectID+=1
			self.removeFileCtg()

	def removeFileCtg(self):
		if self.unplacednum >0:
			if self.remove==1 or self.remove==2:
				if os.path.exists(self.projectData):
					commondline="rm -rf "+self.projectData
					os.system(commondline)
		else:
			if self.remove==1:
				if os.path.exists(self.project):
					commondline="rm -rf "+self.project
					os.system(commondline)

	def skipRound(self):
		logFile=open(self.projectlog,'r')
		for row in logFile:
			if row[0]!='#':
				row1=row.rstrip().split('\t')
				if row1!=[]:
					if row1[0]=='projectID':
						self.projectID=int(row1[1])
					elif row1[0]=='projectOut':
						self.projectOut=row1[1]
					elif row1[0]=='projectData':
						self.projectData=row1[1]
					elif row1[0]=='projectfasta':
						self.projectfasta=row1[1]
						for gseq in SeqIO.parse(self.projectfasta,'fasta'):
							seqout=gseq.seq
					elif row1[0]=='projectagppwd':
						self.projectagppwd=row1[1]
					elif row1[0]=='projectTerminalCtg':
						self.projectagppwd=row1[1]
					elif row1[0]=='TerminalCtgID':
						if len(row1)!=1:
							self.TerminalCtgID=row1[1].split(";")
						else:
							self.TerminalCtgID==[]
					elif row1[0]=='projectused':
						self.projectused=row1[1]
					elif row1[0]=='ctgUsedID':
						if len(row1)==1:
							self.ctgUsedID=[]
						else:
							self.ctgUsedID=row1[1].split(";")
					elif row1[0]=='flaglist':
						self.flaglist=row1[1].split(";")
					elif row1[0]=='roundInputSeqFile':
						self.roundInputSeqFile=row1[1]
					elif row1[0]=='roundInputSeqID':
						self.roundInputSeqID=row1[1]
					elif row1[0]=='roundInputSeqLength':
						self.roundInputSeqLength=int(row1[1])
					elif row1[0]=='placedlist':
						self.placedlist=row1[1].split(";")
					elif row1[0]=='projectagp':
						if len(row1)<=2:
							self.projectagp=row1[1]
					elif row1[0]=='DGUsedCtgList':
						if len(row1)==1:
							self.DGUsedCtgList=[]
						else:
							self.DGUsedCtgList=row1[1].split(";")
					elif row1[0]=='projectusedLine':
						self.projectusedLine='\t'.join(row1[1:])+"\n"
					elif row1[0]=='projectSeq':
						self.projectSeq=row1[1]
					elif row1[0]=='leftTag':
						if len(row1)==1:
							self.leftTag=''
						else:
							self.leftTag=row1[1]
					elif row1[0]=='rightTag':
						if len(row1)==1:
							self.rightTag=''
						else:
							self.rightTag=row1[1]
					elif row1[0]=='unplacedlist':
						if len(row1)==1:
							self.unplacedlist=[]
						else:
							self.unplacedlist=row1[1].split(";")
					elif row1[0]=='unplacednum':
						self.unplacednum=int(row1[1])
		logFile.close()
		DGusetaglist=self.DGUsedCtgList

		file1=open(self.used,'r')
		file2=open(self.projectOut+'/temp','w')
		for row in file1:
			if row[0]!='#':
				row1=row.split('\t')
				if row1[0] in DGusetaglist:
					if self.rightTag=="**" and self.leftTag=="**":
						l=row1[0]+"\t*\t*\t"+row1[-1]
					else:
						l=row1[0]+"\t"+self.rightTag+"\t"+self.leftTag+"\t"+row1[-1]
					file2.writelines(l)
				else:
					file2.writelines(row)
			else:
				file2.writelines(row)
		file1.close()

		file1=open(self.agppwd,'a')
		l="DG"+str(self.projectID)+"\t"+self.projectagp+"\n"
		file1.writelines(l)
		file1.close()

		if self.rightTag=='**' and self.leftTag=='**':
			self.setScaffold(seqout)
			self.scaffoldID+=1
			self.Lastround=True
		else:
			self.Lastround=False
		file2.writelines(self.projectusedLine)
		file2.close()
		file1=open(self.projectSeq,'a')
		for gseq in SeqIO.parse(self.projectfasta,'fasta'):
			l='>'+gseq.description+"\n"+gseq.seq+"\n"
			file1.writelines(l)
		file1.close()
		commondline='mv '+self.projectOut+'/temp '+self.used
		os.system(commondline)

		file1=open(self.unplaced,'r')
		file2=open(self.projectOut+'/temp','w')
		unplacedn=0
		self.unplacedlist=[]
		for row in file1:
			if row.rstrip() not in self.placedlist:
				file2.writelines(row)
				unplacedn+=1
				self.unplacedlist.append(row.rstrip())
		file1.close()
		file2.close()
		commondline='mv '+self.projectOut+'/temp '+self.unplaced
		os.system(commondline)
		self.unplacednum=unplacedn

	def elongationDG(self):
		logFile=open(self.projectlog,'w')
		LogLine="\n\n*****************\n\nprojectID\t"+str(self.projectID)+"\n"
		LogLine+="projectOut\t"+str(self.projectOut)+"\n"
		LogLine+="projectData\t"+str(self.projectData)+"\n"
		LogLine+="projectfasta\t"+str(self.projectfasta)+"\n"
		LogLine+="projectagppwd\t"+str(self.projectagppwd)+"\n"
		print (LogLine)
		logline=self.elongationDGInit()
		LogLine+=logline
		self.setGapFiller()

		logline=self.elongationDGNext()
		LogLine+=logline
		print (LogLine)
		logFile.writelines(LogLine)
		logFile.close()

	def elongationDGNext(self):
		file1=open(self.projectagp,'w')
		ltused=[]
		ltunplaced=[]
		logline='projectagp\t'+self.projectagp+"\n"
		if len(self.outdict)==1:
			for k,v in self.outdict.items():
				filet=open(v.agp,'r')
				rownum=0
				if k=='left':
					s2=-1
					e2=-1
					stag='*'
				else:
					s1=-1
					e1=-1
					etag='*'
				for gseq in SeqIO.parse(v.initialSeq,'fasta'):
					InitialSeq=gseq
				for agpl in filet:
					rownum+=1
					agpl1=agpl.split('\t')
					if 'Right' not in agpl1[5] and 'Left' not in agpl1[5]:
						lt="DG"+str(self.projectID)+"\t"+"\t".join(agpl1[1:])
						ltused.append(agpl1[5])
						ltunplaced.append(agpl1[5])
					else:
						lt="DG"+str(self.projectID)+"\t"+"\t".join(agpl1[1:5])+"\t"+"DG"+str(self.projectID)+"-"+"\t".join(agpl1[5:])
						nametemp=agpl1[5].split("DG"+str(self.projectID))[-1]
						ltused.append(InitialSeq.id+nametemp)
					if k=='left' and InitialSeq.id in agpl:
						s2=int(agpl1[1])-1
						e2=int(agpl1[2])
					if k=='right' and InitialSeq.id in agpl:
						s1=int(agpl1[1])-1
						e1=int(agpl1[2])
					file1.writelines(lt)
					logline=logline+'projectagp\t'+lt
					print (agpl,InitialSeq.id,k=='Left',k, InitialSeq.id in agpl)
			file1.close()
			if rownum==1 and k=='left':
				etag='*'
			else:
				etag=''
			if rownum==1 and k=='right':
				stag='*'
			else:
				stag=''
			GF=v
			for gseq in SeqIO.parse(GF.Elongation.finalSeq,'fasta'):
				GFseq=gseq
				seqoutt=gseq.seq
		else:
			sGF=self.outdict['left']
			eGF=self.outdict['right']
			for gseq in SeqIO.parse(sGF.Elongation.finalSeq,'fasta'):
				sGFseq=gseq
			for gseq in SeqIO.parse(eGF.Elongation.finalSeq,'fasta'):
				eGFseq=gseq
			j=1
			file2=open(eGF.agp,'r')
			s1=-1
			e1=-1
			rownum=0
			lt1=[]
			for gseq in SeqIO.parse(eGF.initialSeq,'fasta'):
				InitialSeq=gseq
			for agpl in file2:
				rownum+=1
				agpl1=agpl.split('\t')
				if 'Right' not in agpl1[5]:
					if InitialSeq.id in agpl:
						lt1=agpl1
						endlen1=int(agpl1[6])-1
					else:
						lt="DG"+str(self.projectID)+"\t"+"\t".join(agpl1[1:])
						logline=logline+'projectagp\t'+lt
						file1.writelines(lt)
						j+=1
					ltused.append(agpl1[5])
					ltunplaced.append(agpl1[5])
				else:
					lt="DG"+str(self.projectID)+"\t"+"\t".join(agpl1[1:5])+"\t"+"DG"+str(self.projectID)+"-"+"\t".join(agpl1[5:])
					nametemp=agpl1[5].split("DG"+str(self.projectID))[-1]
					ltused.append(InitialSeq.id+nametemp)
					file1.writelines(lt)
					j+=1
					logline=logline+'projectagp\t'+lt
				if InitialSeq.id in agpl:
					s1=int(agpl1[1])-1
					e1=int(agpl1[2])
			file2.close()
			if rownum==1:
				stag='*'
			else:
				stag=''
			
			file2=open(sGF.agp,'r')
			if s1==-1:
				endlen=0
			else:
				if endlen1!=0:
					endlen=len(eGFseq.seq[:s1])-endlen1
				else:
					endlen=len(eGFseq.seq[:s1])
			s2=-1
			e2=-1
			rownum=0
			for gseq in SeqIO.parse(eGF.initialSeq,'fasta'):
				InitialSeq=gseq
			for agpl in file2:
				rownum+=1
				agpl1=agpl.split('\t')
				if InitialSeq.id in agpl:
					s2=int(agpl1[1])-1
					e2=int(agpl1[2])
					if lt1!=[]:
						a=int(lt1[7])-int(agpl1[7])
					else:
						a=0
					if s1==-1 and e1==-1:
						lt="DG"+str(self.projectID)+"\t"+str(int(agpl1[1])+endlen-a)+"\t"+str(int(agpl1[2])+endlen-a)+"\t"+str(j)+"\t"+"\t".join(agpl1[4:])
						if 'Left' not in agpl1[5]:
							ltunplaced.append(agpl1[5])
						ltused.append(agpl1[5])
					else:
						lt="DG"+str(self.projectID)+"\t"+lt1[1]+"\t"+str(int(lt1[2])-a)+"\t"+str(j)+"\tw\t"+lt1[5]+"\t"+lt1[6]+"\t"+str(int(lt1[7])-a)+"\t"+lt1[8]
					file1.writelines(lt)
					logline=logline+'projectagp\t'+lt
					j+=1
				else:
					if 'Left' not in agpl1[5]:
						lt="DG"+str(self.projectID)+"\t"+str(int(agpl1[1])+endlen)+"\t"+str(int(agpl1[2])+endlen)+"\t"+str(j)+"\t"+"\t".join(agpl1[4:])
						ltused.append(agpl1[5])
						ltunplaced.append(agpl1[5])
					else:
						lt="DG"+str(self.projectID)+"\t"+str(int(agpl1[1])+endlen)+"\t"+str(int(agpl1[2])+endlen)+"\t"+str(j)+"\t"+agpl1[4]+"\tDG"+str(self.projectID)+"-"+"\t".join(agpl1[5:])
						nametemp=agpl1[5].split("DG"+str(self.projectID))[-1]
						ltused.append(InitialSeq.id+nametemp)
					file1.writelines(lt)
					logline=logline+'projectagp\t'+lt
					j+=1
			file2.close()
			file1.close()
			if rownum==1:
				etag='*'
			else:
				etag=''

		file1=open(self.used,'r')
		file2=open(self.projectOut+'/temp','w')

		if len(ltused)!=0:
			if 'noExtensionContigsorReads' in ltused[0] or 'noNewExtensionReads' in ltused[0] or stag=='*' or 'reachMaximumLength' in ltused[0]:
				stag='*'
			else:
				stag=''
			if 'noExtensionContigsorReads' in ltused[-1] or 'noNewExtensionReads' in ltused[-1] or etag=='*' or 'reachMaximumLength' in ltused[-1]:
				etag='*'
			else:
				etag=''
		if len(self.outdict)!=1:
			if s1==-1 and e1==-1:
				stag='*'
			if s2==-1 and e2==-1:
				etag='*'
		else:
			if 'right' in self.outdict:
				print (s1,e1)
				if s1==-1 and e1==-1:
					stag='*'
				etag='*'
			if 'left' in self.outdict:
				print (s2,e2)
				if s2==-1 and e2==-1:
					etag='*'
				stag='*'
		print (stag,etag,"1,stag,etag!!!!!\n")
		self.placedlist=ltunplaced
		logline=logline+"placedlist\t"+";".join(ltunplaced)+"\n"
		DGusetaglist=[]
		for row in file1:
			print (row)
			row1=row.split('\t')
			if row1[0] in ltused:
				stt=row1[1]
				ett=row1[2]
				if stt=='*':
					stag='*'
				if ett=='*':
					etag='*'
				DGn=row1[0]
				DGusetaglist.append(DGn)
				setDGusetaglist=sorted(list(set(DGusetaglist)))
				while DGn!='':
					file3=open(self.used,'r')
					for r3 in file3:
						r31=r3.rstrip().split('\t')
						if r31[0]==DGn:
							DGn1=DGn
							DGn=''
							r32=r31[3].split(';')
							for dgn in r32:
								print (dgn,'DG' in dgn,r31[0] not in dgn,r31[0])
								if 'DG' in dgn and r31[0] not in dgn:
									DGn=dgn
									if DGn not in DGusetaglist:
										DGusetaglist.append(DGn)
					file3.close()
					if sorted(list(set(DGusetaglist)))==setDGusetaglist:
						DGn=''
					else:
						print (DGusetaglist)
						setDGusetaglist=sorted(list(set(DGusetaglist)))
		print (stag,etag,"2,stag,etag!!!!!\n")
		self.DGUsedCtgList=DGusetaglist
		logline=logline+"DGUsedCtgList\t"+";".join(self.DGUsedCtgList)+"\n"
		file1.close()
		file1=open(self.used,'r')
		file2=open(self.projectOut+'/temp','w')
		for row in file1:
			row1=row.split('\t')
			if row1[0] in DGusetaglist:
				l=row1[0]+"\t"+stag+"\t"+etag+"\t"+row1[-1]
				file2.writelines(l)
			else:
				file2.writelines(row)
		file1.close()
		if stag=='*' and etag=='*':
			stag='**'
			etag='**'
		print (stag,etag,"3,stag,etag!!!!!\n")
		ltusedl="DG"+str(self.projectID)+"\t"+stag+"\t"+etag+"\t"+";".join(ltused)+"\n"
		logline=logline+'projectusedLine\t'+ltusedl+"\n"
		self.projectusedLine=ltusedl
		file3=open(self.projectused,'w')
		file2.writelines(ltusedl)
		file3.writelines(ltusedl)
		file2.close()
		file3.close()
		
		commondline='mv '+self.projectOut+'/temp '+self.used
		os.system(commondline)
		
		file1=open(self.projectSeq,'a')
		if len(self.outdict)==1:
			ChrID=self.roundInputSeq.description.split('\t')[-1]
			l=">DG"+str(self.projectID)+"\t"+";".join(ltused)+"\t"+ChrID+"\n"+GFseq.seq+"\n"
			seqout=GFseq.seq
		else:
			if s2==-1 and e2==-1:
				seqoutt=eGFseq.seq
				ChrID=self.roundInputSeq.description.split('\t')[-1]
			else:
				seqoutt=eGFseq.seq[:e1-a]+sGFseq.seq[e2:]
				ChrID=self.roundInputSeq.description.split('\t')[-1]
			l=">DG"+str(self.projectID)+"\t"+";".join(ltused)+"\t"+ChrID+"\n"+seqoutt+"\n"
		logline=logline+'projectSeq\t'+self.projectSeq+"\n"
		file1.writelines(l)
		file1.close()
		file3=open(self.projectfasta,'w')
		file3.writelines(l)
		file3.close()

		file1=open(self.agppwd,'a')
		l="DG"+str(self.projectID)+"\t"+self.projectagp+"\n"
		file1.writelines(l)
		file1.close()
		file3=open(self.projectagppwd,'w')
		file3.writelines(l)
		file3.close()
		self.leftTag=etag
		self.rightTag=stag
		logline=logline+'leftTag\t'+self.leftTag+"\n"
		logline=logline+'rightTag\t'+self.rightTag+"\n"
		if stag=='**' and etag=='**':
			self.setScaffold(seqoutt)
			self.scaffoldID+=1
			self.Lastround=True
		else:
			self.Lastround=False

		file1=open(self.unplaced,'r')
		file2=open(self.projectOut+'/temp','w')
		unplacedn=0
		self.unplacedlist=[]
		for row in file1:
			row1=row.rstrip().split('\t')
			if row1[0] not in ltunplaced:
				file2.writelines(row)
				unplacedn+=1
				self.unplacedlist.append(row.rstrip())
		file1.close()
		file2.close()
		logline=logline+"unplacedlist\t"+";".join(self.unplacedlist)+"\n"
		print ("unplacedlist\t"+";".join(self.unplacedlist)+"\n")
		commondline='mv '+self.projectOut+'/temp '+self.unplaced
		os.system(commondline)
		self.unplacednum=unplacedn
		logline=logline+"unplacednum\t"+str(self.unplacednum)+"\n"
		return logline

	def setScaffold(self,seqout):
		file1=open(self.scaffoldSeq,'a')
		l='>Scaffold_'+str(self.scaffoldID)+"\t"+"DG"+str(self.projectID)+"\n"+seqout+"\n"
		file1.writelines(l)
		file1.close()
		file1=open(self.agp,'a')
		DGlist=["DG"+str(self.projectID)]
		DGlist+=self.DGUsedCtgList
		for iDG in DGlist:
			file2=open(self.agppwd,'r')
			for agppwd in file2:
				agppwd1=agppwd.rstrip().split('\t')
				if agppwd1[0]==iDG:
					file3=open(agppwd1[1],'r')
					for l in file3:
						lo="Scaffold_"+str(self.scaffoldID)+"("+"DG"+str(self.projectID)+")"+l
						file1.writelines(lo)
					file3.close()
			file2.close()
		file1.close()


	def setGapFiller(self):
		outdict={}
		for flag1 in self.flaglist:
			for gseq in SeqIO.parse(self.roundInputSeqFile,'fasta'):
				leftID=''
				rightID=''
				if 'DG' not in gseq.id:
					chrID=gseq.description.split('\t')[-1]
					contigID=int(gseq.id.split('-')[-1])
					for gseq1 in SeqIO.parse(self.projectTerminalCtg,'fasta'):
						chrID1=gseq1.description.split('\t')[-1]
						contigID1=int(gseq1.id.split('-')[-1])
						if chrID1==chrID:
							if contigID1==(contigID+1):
								leftID=gseq1.id
							if contigID1==(contigID-1):
								rightID=gseq1.id
				else:
					chrID=gseq.description.split('\t')[-2]
					print (gseq.description)
					contigID=gseq.description.split('\t')[-3].split(';')
					contigIDright=''
					contigIDleft=''
					if 'Right' not in contigID[0] and len(contigID[0].split('-'))!=1:
						print (contigID[0])
						contigIDright=int(contigID[0].split('-')[1])
					if 'Left' not in contigID[-1] and len(contigID[-1].split('-'))!=1:
						contigIDleft=int(contigID[-1].split('-')[1])
					for gseq1 in SeqIO.parse(self.projectTerminalCtg,'fasta'):
						chrID1=gseq1.description.split('\t')[-1]
						contigID1=int(gseq1.id.split('-')[-1])
						if chrID1==chrID:
							if contigIDleft!='':
								if contigID1==(contigIDleft+1):
									leftID=gseq1.id
							if contigIDright!='':
								if contigID1==(contigIDright-1):
									rightID=gseq1.id
					
				
					
			if flag1=='left':
				projectoutEX=self.projectOut+"/DG"+str(self.projectID)+"-Left"
			else:
				projectoutEX=self.projectOut+"/DG"+str(self.projectID)+"-Right"
			if not os.path.exists(projectoutEX):
				os.makedirs(projectoutEX)
			
			self.ctgsEdgeFile=self.projectData+"/ctgsEdge.fa"
			ft=open(self.ctgsEdgeFile,'w')
			seedlen=int(self.seedLen)
			if leftID!='' and flag1=='left':
				print (flag1,leftID)
				for gseq in SeqIO.parse(self.projectTerminalCtg,'fasta'):
					if gseq.id==leftID:
						if len(gseq.seq)<=seedlen:
							l='>'+gseq.id+"-edge-forward_com\n"+gseq.seq+"\n"
							ft.writelines(l)
						else:
							l='>'+gseq.id+"-edge-forward\n"+gseq.seq[:seedlen]+"\n"
							ft.writelines(l)
				ft.close()
			elif rightID!='' and flag1!='left':
				print (flag1,rightID)
				for gseq in SeqIO.parse(self.projectTerminalCtg,'fasta'):
					if gseq.id==rightID:
						if len(gseq.seq)<=seedlen:
							l='>'+gseq.id+"-edge-forward_com\n"+gseq.seq+"\n"
							ft.writelines(l)
						else:
							l='>'+gseq.id+"-edge-forward\n"+gseq.seq[-seedlen:]+"\n"
							ft.writelines(l)
				ft.close()
			else:
				for gseq in SeqIO.parse(self.projectTerminalCtg,'fasta'):
					if len(gseq.seq)<=seedlen:
						l='>'+gseq.id+"-edge-forward_com\n"+gseq.seq+"\n"
						ft.writelines(l)
						l='>'+gseq.id+"-edge-reverse_com\n"+gseq.seq.reverse_complement()+"\n"
						ft.writelines(l)
						
					else:
						if flag1=='left':
							l='>'+gseq.id+"-edge-forward\n"+gseq.seq[:seedlen]+"\n"
							ft.writelines(l)
							l='>'+gseq.id+"-edge-reverse\n"+gseq.seq.reverse_complement()[:seedlen]+"\n"
							ft.writelines(l)
						
						else:
							l='>'+gseq.id+"-edge-forward\n"+gseq.seq[-seedlen:]+"\n"
							ft.writelines(l)
							l='>'+gseq.id+"-edge-reverse\n"+gseq.seq.reverse_complement()[-seedlen:]+"\n"
							ft.writelines(l)
				ft.close()	
				
			if flag1=='left':
				gfout=GapFiller([self.mode,self.remove,self.thread,self.reads,projectoutEX,self.roundInputSeqFile,self.ctgsEdgeFile,flag1,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict,self.maxReadsLen,self.seedLen])
			else:
				gfout=GapFiller([self.mode,self.remove,self.thread,self.reads,projectoutEX,self.ctgsEdgeFile,self.roundInputSeqFile,flag1,self.edge,self.filterDepth,self.MaximunExtensionLength,self.readsDict,self.maxReadsLen,self.seedLen])
			outdict[flag1]=gfout
			if gfout.Elongation.roundResult.ExtensionReads.note=='' and 'No extension contigs or reads found' not in gfout.Elongation.roundResult.ExtensionContigs.selectContigNote and "Reach the maximum Length" not in gfout.Elongation.roundResult.ExtensionContigs.selectContigNote:
				seedlen=int(self.seedLen)
				if gfout.Elongation.roundResult.roundOutput.linkedSequenceNote!='':
					for gseq in SeqIO.parse(gfout.Elongation.roundResult.roundOutput.linkedSequence,'fasta'):
						l0=gseq.description
						l1=l0.split('\t')
						for l11 in l1:
							l2=l11.split(':')
							if l2[0]=="Aln":
								closectg=l2[1].split(';')
								TS=closectg[-2]
					TS1=TS.split("-edge-")
					for gseq in SeqIO.parse(self.ctgSeq,'fasta'):
						if gseq.id==TS1[0]:
							TSSeq=gseq
						if gseq.id==self.roundInputSeq.id:
							ISSeq=gseq
					for gseq in SeqIO.parse(gfout.Elongation.finalSeq,'fasta'):
						GSseq=gseq

					if 'reverse' in TS1[1]:
						TSSeqS=TSSeq.seq.reverse_complement()
					else:
						TSSeqS=TSSeq.seq
					if flag1=='left':
						if 'com' not in TS1[1]:
							SeqFinal=GSseq.seq+TSSeqS[int(seedlen):]
							addlen=len(TSSeqS[int(seedlen):])
						else:
							SeqFinal=GSseq.seq
							addlen=0
					else:
						if 'com' not in TS1[1]:
							SeqFinal=TSSeqS[:-int(seedlen)]+GSseq.seq
							addlen=len(TSSeqS[:-int(seedlen)])
						else:
							SeqFinal=GSseq.seq
							addlen=0
	
					filet=open(gfout.agp,'r')
					lage=''
					for row in filet:
						lt=''
						row1=row.rstrip().split('\t')
						if flag1=='left':
							if TS1[0] in row:
								lt=row1[0]+"\t"+row1[1]+"\t"+str(len(SeqFinal))+"\t"+row1[3]+"\t"+row1[4]+"\t"+TS1[0]+"\t"+row1[6]+"\t"+str(int(row1[7])+addlen)
								if 'reverse' in TS1[1]:
									lt=lt+"\t-\n"
								else:
									lt=lt+"\t+\n"
								lage+=lt
							else:
								lage+=row
						else:
							if TS1[0] in row:
								lt=row1[0]+"\t"+row1[1]+"\t"+str(int(row1[7])+addlen)+"\t"+row1[3]+"\t"+row1[4]+"\t"+TS1[0]+"\t"+row1[6]+"\t"+str(int(row1[7])+addlen)
								if 'reverse' in TS1[1]:
									lt=lt+"\t-\n"
								else:
									lt=lt+"\t+\n"
								lage+=lt
							else:
								if row1[5]==GSseq.id:
									lt=row1[0]+"\t"+str(int(row1[1])+addlen)+"\t"+str(int(row1[2])+addlen)+"\t"+row1[3]+"\t"+row1[4]+"\t"+row1[5]+"\t"+str(int(row1[6])+addlen)+"\t"+str(int(row1[7])+addlen)+"\t+\n"
								else:
									lt=row1[0]+"\t"+str(int(row1[1])+addlen)+"\t"+str(int(row1[2])+addlen)+"\t"+"\t".join(row1[3:7])+"\t"+str(len(self.roundInputSeq.seq))+"\t+\n"
								lage+=lt
						
						if len(SeqFinal)!=len(GSseq.seq)+addlen:
							print ('linkage err')
							sys.exit()
					filet.close()
					filet=open(gfout.agp,'w')
					filet.writelines(lage)
					filet.close()
					filet=open(gfout.Elongation.finalSeq,'w')
					lo='>'+GSseq.description+"\n"+SeqFinal+"\n"
					filet.writelines(lo)
					filet.close()
			self.outdict=outdict

	def elongationDGInit(self):
		self.elongationDGInitLog=self.projectOut+"/DG"+str(self.projectID)+".input.log"
		if os.path.exists(self.projectTerminalCtg)!=True or os.path.getsize(self.projectTerminalCtg)==0 or os.path.getsize(self.elongationDGInitLog)==0 or os.path.exists(self.elongationDGInitLog)!=True:
			file1=open(self.projectTerminalCtg,'w')
			fileDGn=open(self.used,'r')
			logfile=open(self.elongationDGInitLog,'w')
			self.ctgUsedID=[]
			self.TerminalCtgID=[]
			for i in fileDGn:
				if i[0]=="#":
					continue
				i1=i.rstrip().split('\t')
				print (i1)
				ctgn=i1[3].split(';')
				print (ctgn)
				for n1 in ctgn:
					if 'DG' not in n1:
						if 'Left' not in n1 and 'Right' not in n1:
							self.ctgUsedID.append(n1)
					else:
						if n1!=i1[0]:
							self.ctgUsedID.append(n1)
				if '*' in i1[1] and '*' in i1[2]:
					continue
				else:
					for gseq in SeqIO.parse(self.projectSeq,'fasta'):
						if gseq.id==i1[0]:
							l='>'+gseq.description+"\t"+i1[1]+";"+i1[2]+"\n"+gseq.seq+"\n"
							file1.writelines(l)

			fileDGn.close()
			print (self.ctgUsedID,'self.ctgUsedID')
			for gseq in SeqIO.parse(self.ctgSeq,'fasta'):
				if gseq.id not in self.ctgUsedID:
					l='>'+gseq.description+"\n"+gseq.seq+"\n"
					file1.writelines(l)
			file1.close()
			self.FindInputSeq()

			LogLine="projectTerminalCtg\t"+str(self.projectTerminalCtg)+"\n"
			LogLine+="TerminalCtgID\t"+";".join(self.TerminalCtgID)+"\n"
			LogLine+="projectused\t"+str(self.projectused)+"\n"
			LogLine+="ctgUsedID\t"+";".join(self.ctgUsedID)+"\n"
			LogLine+="flaglist\t"+";".join(self.flaglist)+"\n"

			LogLine+="roundInputSeqFile\t"+str(self.roundInputSeqFile)+"\n"
			LogLine+="roundInputSeqID\t"+str(self.roundInputSeq.id)+"\n"
			LogLine+="roundInputSeqLength\t"+str(len(self.roundInputSeq.seq))+"\n"
			
			logfile.writelines(LogLine)
			logfile.close()
			return LogLine
		else:
			logfile=open(self.elongationDGInitLog,'r')
			LogLine=''
			for row in logfile:
				LogLine+=row
				row1=row.rstrip().split('\t')
				if row1[0]=='projectTerminalCtg':
					self.projectTerminalCtg=row1[1]
				elif row1[0]=='TerminalCtgID':
					if len(row1)>1:
						self.TerminalCtgID=row1[1].split(';')
					else:
						self.TerminalCtgID=[]
				elif row1[0]=='projectused':
					self.projectused=row1[1]
				elif row1[0]=='ctgUsedID':
					if len(row1)>1:
						self.ctgUsedID=row1[1].split(';')
					else:
						self.ctgUsedID=[]
				elif row1[0]=='flaglist':
					self.flaglist=row1[1].split(';')
				elif row1[0]=='roundInputSeqFile':
					self.roundInputSeqFile=row1[1]
					for gseq in SeqIO.parse(self.roundInputSeqFile,'fasta'):
						self.roundInputSeq=gseq
			logfile.close()
			return LogLine
					
	def FindInputSeq(self):
		self.roundInputSeqFile=self.projectData+"/DG"+str(self.projectID)+".input.fa"
		lenm=0
		for gseq in SeqIO.parse(self.projectTerminalCtg,'fasta'):
			if len(gseq.seq)>=lenm:
				if len(gseq.seq)==lenm:
					ID1=int(filter(str.isdigit, gseq.id))
					ID2=int(filter(str.isdigit, self.roundInputSeq.id))
					if ID1<ID2:
						lenm=len(gseq.seq)
						self.roundInputSeq=gseq
				else:
				
					lenm=len(gseq.seq)
					self.roundInputSeq=gseq
		ft=open(self.roundInputSeqFile,'w')
		seedlen=self.seedLen
		l='>'+self.roundInputSeq.description+"\n"+self.roundInputSeq.seq+"\n"
		ft.writelines(l)
		ft.close()
		des=self.roundInputSeq.description
		print (self.roundInputSeq.description,'des')
		des1=des.rstrip().split('\t')
		ft=open(self.used,'r')
		rtagt=''
		ltagt=''
		for row in ft:
			row1=row.rstrip().split('\t')
			if row1[-1]==des1[1]:
				rtagt=row1[1]
				ltagt=row1[2]
				print (row1)
				print (rtagt)
				print (ltagt)
		ft.close()
		if len(des1)==1:
			self.flaglist=['left','right']
		else:
			self.flaglist=[]
			des2=des1[1].split(';')
			if 'noExtensionContigsorReads' not in des2[-1] and 'noNewExtensionReads' not in des2[-1] and ltagt=='' and 'reachMaximumLength' not in des2[-1]:
				self.flaglist.append('left')
			if 'noExtensionContigsorReads' not in des2[0] and 'noNewExtensionReads' not in des2[0] and rtagt=='' and 'reachMaximumLength' not in des2[0]:
				self.flaglist.append('right')
		print (self.flaglist)
		ft=open(self.projectData+'/temp','w')
		for gseq in SeqIO.parse(self.projectTerminalCtg,'fasta'):
			if gseq.id!=self.roundInputSeq.id and gseq.id not in self.ctgUsedID:
				l='>'+gseq.description+"\n"+gseq.seq+"\n"
				self.TerminalCtgID.append(gseq.id)
				ft.writelines(l)
		ft.close()
		commondline='mv '+self.projectData+'/temp '+self.projectTerminalCtg
		os.system(commondline)

	def roundInit(self):
		self.projectOut=self.project+"/DG"+str(self.projectID)
		if not os.path.exists(self.projectOut):
			os.makedirs(self.projectOut)
		self.projectData=self.projectOut+"/data"
		if not os.path.exists(self.projectData):
			os.makedirs(self.projectData)
		self.projectTerminalCtg=self.projectOut+"/data/terminalCtg.fa"
		self.projectused=self.projectOut+"/DG"+str(self.projectID)+".used.txt"
		self.projectfasta=self.projectOut+"/DG"+str(self.projectID)+".fa"
		self.projectagppwd=self.projectOut+"/DG"+str(self.projectID)+".agp.path.txt"
		self.projectlog=self.projectOut+"/DG"+str(self.projectID)+".log"
		self.projectagp=self.projectOut+"/DG"+str(self.projectID)+".agp"
		
		
