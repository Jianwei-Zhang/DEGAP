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
import GapFillerClass

class FindExtensionContigs(object):
	def __init__(self,ExtensionReads):
		self.extensionReads=ExtensionReads
		self.log=self.extensionReads.roundInput.elongation.roundDir+"/extensionConthgs.log"
		self.extensionSequence=self.extensionReads.roundInput.elongation.roundDir+"/extensionSequence."+self.extensionReads.roundInput.elongation.base.tag+".fa"
		self.extensionContigID=[]
		self.extensionLength=0
		if os.path.exists(self.extensionSequence)==True and os.path.getsize(self.extensionSequence)!=0:
			self.readlog()
		else:
			logfilet=open(self.log,'w')

			self.hifiasm()
			logLine='hifiasmCommand\t'+self.hifiasmCommand+"\n"
			logLine+='hifiasmDir\t'+self.hifiasmDir+"\n"
			logLine+='hifiasmResult\t'+self.hifiasmResult+"\n"
			listt=[]
			for k,v in self.hifiasmResultDict.items():
				l=k+";"+';'.join(v)
				listt.append(l)
			logLine+='hifiasmResultDict\t'+':'.join(listt)+"\n"

			if os.path.exists(self.hifiasmResult)==True and os.path.getsize(self.hifiasmResult)!=0:
				self.extensionSeqNote='hifiasm'
				self.readCommonResult='None'
				self.extensionContigs=self.hifiasmResult
			else:
				self.readCommon()
				self.extensionContigs=self.readCommonResult
			logLine+='extensionSeqNote\t'+self.extensionSeqNote+"\n"
			logLine+='readCommonResult\t'+self.readCommonResult+"\n"
			logLine+='extensionContigs\t'+self.extensionContigs+"\n"
			self.selectContigNote=''
			self.selectContigNote=self.selectContig()
			if 'hifiasmHasNoAccurateAssembly' in self.selectContigNote:
				print (self.selectContigNote)
				self.readCommon()
				self.extensionContigs=self.readCommonResult
				self.selectContigNote=self.selectContig()
				if len(self.selectExtensionContigsAln)==0:
					self.selectContigNote+='No extension contigs or reads found\n'
				logLine+='selectContigNote\t'+self.selectContigNote+"\n"

			logLine+='potentialExtensionContigsAln\t'+self.potentialExtensionContigsAln+"\n"
			logLine+='selectContigIdentity\t'+str(self.selectContigIdentity)+"\n"
			logLine+='selectContigDistance\t'+str(self.selectContigDistance)+"\n"
			logLine+='selectContigAlnLength\t'+str(self.selectContigAlnLength)+"\n"
			logLine+='contigAlnMerge\t'+str(self.contigAlnMerge)+"\n"
			logLine+='contigAlnMergeIdentity\t'+str(self.contigAlnMergeIdentity)+"\n"
			logtmp=''
			for k,v in self.potentialExtensionContigsAlnDict.items():
				listRefForward,listQueryForward,listRefReverse,listQueryReverse,qlen=v
				l=k+"-"
				for i in listRefForward:
					s,e=i
					lt=str(s)+","+str(e)+";"
					l+=lt
				l+="-"
				for i in listQueryForward:
					s,e=i
					lt=str(s)+","+str(e)+";"
					l+=lt
				l+="-"
				for i in listRefReverse:
					s,e=i
					lt=str(s)+","+str(e)+";"
					l+=lt
				l+='-'
				for i in listQueryReverse:
					s,e=i
					lt=str(s)+","+str(e)+";"
					l+=lt
				l+='-'+str(qlen)+"\t"
				logtmp+=l
			logLine+='potentialExtensionContigsAlnDict\t'+logtmp+"\n"

			logtmp=''
			for i in self.selectExtensionContigsAln:
				r,el=i
				l=r.rstrip()+"-"+str(el)+";"
				logtmp+=l
			logLine+='selectExtensionContigsAln\t'+logtmp+"\n"
			logLine+='maxextensionLen\t'+str(self.maxextensionLen)+"\n"
			logLine+='selectContigNote\t'+self.selectContigNote+"\n"
			logtmp=''
			for i in self.FuzzyAln:
				l=i[0].rstrip()+'-'+str(i[1])+";"
				logtmp+=l
			logLine+='FuzzyAln\t'+logtmp+"\n"
			logLine+='extensionContigID\t'+";".join(self.extensionContigID)+"\n"
			logLine+='extensionLength\t'+str(self.extensionLength)+"\n"
			logLine+='extensionSequence\t'+self.extensionSequence+"\n"
			logfilet.writelines(logLine)
			logfilet.close()

	def selectContig(self):
		self.potentialExtensionContigsAln=self.extensionReads.roundInput.elongation.roundDir+"/potentialExtensionContigs.asm."+self.extensionReads.roundInput.elongation.base.tag+"seq.mummer"
		self.potentialExtensionContigsAln=GapFillerClass.mummer(self.extensionReads.roundInput.inputSeq,self.extensionContigs,self.potentialExtensionContigsAln)
		self.selectContigIdentity=99.0
		self.selectContigDistance=10
		self.selectContigAlnLength=1000
		self.contigAlnMerge=False
		self.contigAlnMergeIdentity=100.0

		self.contigAlnFilter(self.selectContigIdentity)

		while len(self.selectExtensionContigsAln)<=0:
			if self.selectContigAlnLength<=500:
				if self.selectContigIdentity<95.0:
					if self.selectContigDistance>self.extensionReads.roundInput.elongation.base.edge or self.selectContigDistance>self.maxextensionLen:
						self.contigAlnMerge=True		
					else:
						self.selectContigDistance+=10
				else:
					self.selectContigIdentity-=1.0
			else:
				self.selectContigAlnLength-=100
			if self.contigAlnMerge==False:
				self.contigAlnFilter(self.selectContigIdentity)
			else:
				break

		if self.contigAlnMerge==True:
			self.selectContigNote+='FuzzyAlignments\n'
			self.selectContigDistance=10
			while len(self.selectExtensionContigsAln)<=0 and self.contigAlnMergeIdentity>=80.0 and self.selectContigDistance<=self.extensionReads.roundInput.elongation.base.edge:
				self.potentialExtensionContigsAln=self.extensionReads.roundInput.elongation.roundDir+"/potentialExtensionContigs.asm."+self.extensionReads.roundInput.elongation.base.tag+"seq.mummer.delta.filter.coords"
				if self.contigAlnMergeIdentity<=80.0:
					if self.selectContigDistance<=self.extensionReads.roundInput.elongation.base.edge:
						self.selectContigDistance+=10
				else:
					self.contigAlnMergeIdentity-=1.0
					print (self.contigAlnMergeIdentity)
				self.contigAlnFilter(self.contigAlnMergeIdentity)
				self.FuzzyAln=self.mergeContig()
				print (self.FuzzyAln)
				self.selectExtensionContigsAln=self.FuzzyAln
				
				#self.contigAlnFilter(self.contigAlnMergeIdentity)
				
			if len(self.selectExtensionContigsAln)<=0 and self.contigAlnMergeIdentity<=80.0:
				self.selectContigNote+='hifiasmHasNoAccurateAssembly Use readsCommonSequence\n'
				return self.selectContigNote
		else:
			self.selectContigNote+='SolidAlignments\n'
			self.FuzzyAln=[]

		if len(self.selectExtensionContigsAln)>1:
			asm=[]
			self.selectContigNote+="\t\tFind maximum identity alignment"
			for i in self.selectExtensionContigsAln:
				aln,mlen=i
				if asm==[]:
					asm.append(i)
				else:
					ia=asm[0][0].rstrip().split('\t')[6]
					i1=aln.rstrip().split('\t')[6]
					if float(i1)>float(ia):
						asm=[]
						asm.append(i)
					elif float(i1)==float(ia):
						asm.append(i)
			if asm[0][0].rstrip().split('\t')[6]==0 or len(asm)!=1:
				self.selectContigNote+="\t\tFind maximum length alignment\n"
				asm1=[]
				for i in asm:
					aln,mlen=i
					if asm1==[]:
						asm1.append(i)
					else:
						if asm1[0][1]<mlen:
							asm1=[]
							asm1.append(i)
			else:
				asm1=asm
			self.selectExtensionContigsAln=asm1
		file1=open(self.extensionSequence,'w')
		for i in self.selectExtensionContigsAln:
			aln,mlen=i
			outseq,extensionseq,noseq,note1,name1,name2=self.extensionSeq(aln)
			l='>ExtensionSequence-'+self.extensionReads.roundInput.elongation.base.name+"\t"+self.extensionReads.roundInput.inputSeedSequence.id+"\t"+name1+"\n"+outseq+"\n"
			file1.writelines(l)
			self.selectContigNote+=note1
			self.extensionContigID=self.extensionContigID+name2
			self.extensionLength=len(outseq)-len(self.extensionReads.roundInput.inputSeedSequence.seq)
		file1.close()
		return self.selectContigNote

	def extensionSeq(self,row):
		seqtseq=self.extensionReads.roundInput.inputSeedSequence
		row1=row.rstrip().split('\t')
		a,b,c,d=int(row1[0]),int(row1[1]),int(row1[2]),int(row1[3])
		e,f=int(row1[4]),int(row1[5])
		name1=row1[-1]
		name2=name1.split('-ovl-')
		list0=[]
		flag=self.extensionReads.roundInput.elongation.base.flag
		note1="selectContigNote\tExtension: "+name1+"\n"
		for gseq in SeqIO.parse(self.extensionContigs,'fasta'):
			if gseq.id==name1:
				if flag!='left':
					if d>=int(row1[8])-10:
						extensionseq=gseq.seq[:c-1]
						note1=note1+'selectContigNote\tExtensionseq alignment:'+row
					else:
						extensionseq1=gseq.seq.reverse_complement()
						nend=len(gseq.seq)-c
						note1=note1+'selectContigNote\tExtensionseq alignment reverse:'+row
						extensionseq=extensionseq1[:nend]
				else:
					if c<=10:
						extensionseq=gseq.seq[d:]
						note1=note1+'selectContigNote\tExtensionseq alignment:'+row
					else:
						extensionseq1=gseq.seq.reverse_complement()
						nstart=len(gseq.seq)-d+1
						extensionseq=extensionseq1[nstart:]
						note1=note1+'selectContigNote\tExtensionseq alignment reverse:'+row
		if flag!='left':
			noseq=seqtseq.seq[int(row1[0])-1:]
			note1=note1+'selectContigNote\t\tSeedsequence length: '+str(len(seqtseq.seq))+"\nselectContigNote\t\tExtensionlength: "+str(len(extensionseq))+"\n"
			outseq=extensionseq+noseq
		else:
			noseq=seqtseq.seq[:int(row1[1])]
			note1=note1+'selectContigNote\t\tSeedsequence length: '+str(len(seqtseq.seq))+"\nselectContigNote\t\tExtensionlength: "+str(len(extensionseq))+"\n"
			outseq=noseq+extensionseq
		note1=note1+'selectContigNote\t\tTotalLength: '+str(len(outseq))
		return outseq,extensionseq,noseq,note1,name1,name2

	def merge(self,list1):
		list2=sorted(list1,key=lambda x:x[0])
		j=0
		while j<len(list2):
			j1=j+1
			while j1<len(list2):
				a,b=list2[j]
				c,d=list2[j1]
				x=range(a,b+1,1)
				y=range(c,d+1,1)
				t1=sorted(list(set(x)&set(y)))
				if len(t1)!=0:
					t2=sorted(list(set(x)|set(y)))
					list2[j]=[t2[0],t2[-1]]
					list2.pop(j1)
					j1-=1
				j1+=1
			j+=1
		list3=[]
		listLen=0
		for i in list2:
			s,e=i
			if s!=e:
				list3.append(i)
				listLen=listLen+(e-s+1)
		return list3,listLen

	def mergeContig(self):
		FuzzyAln=[]
		for k,v in self.potentialExtensionContigsAlnDict.items():
			listRefForward,listQueryForward,listRefReverse,listQueryReverse,QueryLen=v
			refAalnmergeForward,refAlnmergeForwardLen=self.merge(listRefForward)
			refAalnmergeReverse,refAlnmergeReverseLen=self.merge(listRefReverse)
			if refAlnmergeForwardLen>refAlnmergeReverseLen:
				refaln=listRefForward
				queryaln=listQueryForward
				tag='f'
			else:
				refaln=listRefReverse
				queryaln=listQueryReverse
				tag='r'
			refalnmerge,refalnlen=self.merge(refaln)
			if len(refalnmerge)!=1:
				effectiveQueryAln=[]
				effectiveRefAln=[]
				if self.extensionReads.roundInput.elongation.base.flag!='left':
					effectiveRefEdgeAln=refalnmerge[0]
				else:
					effectiveRefEdgeAln=refalnmerge[-1]
				j=0
				s,e=effectiveRefEdgeAln
				while j<len(refaln):
					s1,e1=refaln[j]
					if s1>=s and e1<=e:
						effectiveQueryAln.append(queryaln[j])
					j+=1
				effectiveRefAln=[effectiveRefEdgeAln]
				queryaln=effectiveQueryAln
			else:
				effectiveRefAln=refalnmerge
			queryalnmerge,queryalnlen=self.merge(queryaln)
			effectiveQueryAln=queryalnmerge
			
			refFuzzyAln=[]
			for i in effectiveRefAln:
				s1,e1=i
				if self.extensionReads.roundInput.elongation.base.flag!='left':
					if s1==1:
						refFuzzyAln.append(i)
				else:
					if e1==len(self.extensionReads.roundInput.inputSeedSequence.seq):
						refFuzzyAln.append(i)
			queryFuzzyAln=[]
			for i in effectiveQueryAln:
				s1,e1=i
				if self.extensionReads.roundInput.elongation.base.flag!='left':
					if e1==QueryLen:
						queryFuzzyAln.append(i)
				else:
					if s1==1:
						queryFuzzyAln.append(i)
			if refFuzzyAln!=[] and queryFuzzyAln!=[]:
				for i in refFuzzyAln:
					am,bm=i
					for i1 in queryFuzzyAln:
						cm,dm=i1
						em=bm-am
						if tag=='f':
							fm=dm-cm+1
							lit=[str(am),str(bm),str(cm),str(dm),str(em),str(fm),str(self.contigAlnMergeIdentity),str(0),str(QueryLen),str(0),str(0),str(0),str(k)]
							lit='\t'.join(lit)+'\n'
						else:
							fm=dm-cm+1
							lit=[str(am),str(bm),str(dm),str(cm),str(em),str(fm),str(self.contigAlnMergeIdentity),str(0),str(QueryLen),str(0),str(0),str(0),str(k)]
							lit='\t'.join(lit)+'\n'
						FuzzyAln.append([lit,max(em,fm)])
			else:
				am=effectiveRefAln[0][0]
				bm=effectiveRefAln[-1][1]
				cm=effectiveQueryAln[0][0]
				dm=effectiveQueryAln[-1][1]
				em=bm-am
				if tag=='f':
					fm=dm-cm+1
					lit=[str(am),str(bm),str(cm),str(dm),str(em),str(fm),str(self.contigAlnMergeIdentity),str(0),str(QueryLen),str(0),str(0),str(0),str(k)]
					lit='\t'.join(lit)+'\n'
				else:
					fm=dm-cm+1
					lit=[str(am),str(bm),str(dm),str(cm),str(em),str(fm),str(self.contigAlnMergeIdentity),str(0),str(QueryLen),str(0),str(0),str(0),str(k)]
					lit='\t'.join(lit)+'\n'
				if fm < QueryLen:
					if refFuzzyAln!=[] or queryFuzzyAln!=[]:
						FuzzyAln.append([lit,max(em,fm)])
		self.potentialExtensionContigsAln=self.extensionReads.roundInput.elongation.roundDir+"/potentialExtensionContigs.asm.FuzzyAln.delta.filter.coords"
		file1=open(self.potentialExtensionContigsAln,'w')
		for i in FuzzyAln:
			i1,i2=i
			file1.writelines(i1)
		file1.close()
		return FuzzyAln

	def contigAlnFilter(self,selectIdentity):
		file1=open(self.potentialExtensionContigsAln,'r')
		self.potentialExtensionContigsAlnDict={}
		self.selectExtensionContigsAln=[]
		self.maxextensionLen=0
		for row in file1:
			row1=row.rstrip().split('\t')
			s1,e1,s2,e2=int(row1[0]),int(row1[1]),int(row1[2]),int(row1[3])
			if self.extensionReads.roundInput.elongation.base.flag!='left':
				refDistance=s1-1
				if e2>=s2:
					extensionLen=s2-1
					queryDistance=int(row1[8])-e2
				else:
					extensionLen=int(row1[8])-s2
					queryDistance=e2-1
			else:
				refDistance=int(row1[7])-e1
				if s2<e2:
					extensionLen=int(row1[8])-e2
					queryDistance=s2-1
				else:
					extensionLen=e2-1
					queryDistance=int(row1[8])-s2

			if float(row1[6])>=selectIdentity:
				if row1[12] in self.potentialExtensionContigsAlnDict:
					listRefForward,listQueryForward,listRefReverse,listQueryReverse,Querylen=self.potentialExtensionContigsAlnDict[row1[12]]
				else:
					listRefForward=[]
					listQueryForward=[]
					listRefReverse=[]
					listQueryReverse=[]
				if e2>=s2:
					listRefForward.append([s1,e1])
					listQueryForward.append([s2,e2])
				else:
					listRefReverse.append([s1,e1])
					listQueryReverse.append([e2,s2])
				self.potentialExtensionContigsAlnDict[row1[12]]=[listRefForward,listQueryForward,listRefReverse,listQueryReverse,int(row1[8])]
				if extensionLen>self.maxextensionLen:
					self.maxextensionLen=extensionLen
				if queryDistance<=self.selectContigDistance and refDistance<=self.selectContigDistance:
					maxAlnLen=max([int(row1[4]),int(row1[5])])
					if maxAlnLen>=self.selectContigAlnLength and extensionLen>=0:
						self.selectExtensionContigsAln.append([row,extensionLen])
		file1.close()

	def readCommon(self):
		self.readCommonResult=self.extensionReads.roundInput.elongation.roundDir+"/readsCommonSequence.fasta"
		mummerout=GapFillerClass.mummer(self.extensionReads.roundInput.inputSeq,self.extensionReads.extensionReads,self.extensionReads.roundInput.elongation.roundDir+'/temp3')
		ft3=open(mummerout,'r')
		readlisttemp=[]
		for r in ft3:
			r1=r.rstrip().split('\t')
			if r1[-1] not in readlisttemp:
				readlisttemp.append(r1[-1])
		ft3.close()
		file1=open(self.readCommonResult,'w')
		if self.extensionReads.extensionReadsNum >1:
			self.extensionSeqNote='readsCommon'
			listCommon=[]
			j=1
			for gseq1 in SeqIO.parse(self.extensionReads.extensionReads,'fasta'):
				if gseq1.id not in readlisttemp:
					continue
				ft1=open(self.extensionReads.roundInput.elongation.roundDir+'/temp1','w')
				l='>'+gseq1.id+'\n'+gseq1.seq+"\n"
				ft1.writelines(l)
				ft1.close()
				readSeq1=gseq1.seq
				listCommon.append(gseq1.id)

				for gseq2 in SeqIO.parse(self.extensionReads.extensionReads,'fasta'):
					if gseq2.id not in readlisttemp:
						continue
					if gseq1.id != gseq2.id:
						if gseq2.id in listCommon:
							continue
						else:
							j=j+1
							ft2=open(self.extensionReads.roundInput.elongation.roundDir+'/temp2','w')
							readSeq2=gseq2.seq
							l='>'+gseq2.id+'\n'+gseq2.seq+"\n"
							ft2.writelines(l)
							ft2.close()

							mummerout=GapFillerClass.mummer(self.extensionReads.roundInput.elongation.roundDir+'/temp2',self.extensionReads.roundInput.elongation.roundDir+'/temp1',self.extensionReads.roundInput.elongation.roundDir+'/temp3')
							ft3=open(mummerout,'r')
							maxCommonLen=0
							mummeresult=[]
							for r in ft3:
								r1=r.rstrip().split('\t')
								if int(r1[0])==1 or int(r1[0])==int(r1[7]) or int(r1[1])==1 or int(r1[1])==int(r1[7]):
							
									if int(r1[2])==1 or int(r1[2])==int(r1[8]) or int(r1[3])==1 or int(r1[3])==int(r1[8]):
										maxAlnlen=max(float(r1[5]),float(r1[4]))
										if maxAlnlen>=maxCommonLen:
											maxCommonLen=maxAlnlen
											mummeresult=r1
							ft3.close()

							if mummeresult!=[]:
								if int(mummeresult[4])>int(mummeresult[5]):
									seq1=readSeq2[int(mummeresult[0])-1:int(mummeresult[1])]
								else:
									if int(mummeresult[2])<int(mummeresult[3]):
										seq1=readSeq1[int(mummeresult[2])-1:int(mummeresult[3])]
									else:
										seq1=readSeq1[int(mummeresult[3])-1:int(mummeresult[2])]
								l='>'+gseq1.id+"-ovl-"+gseq2.id+"\tcommon-"+str(len(seq1))+'\t'+r+seq1+"\n"
								file1.writelines(l)
							commendline='rm '+self.extensionReads.roundInput.elongation.roundDir+'/temp2'
							os.system(commendline)
							commendline='rm '+self.extensionReads.roundInput.elongation.roundDir+'/temp3*'
							os.system(commendline)
				commendline='rm '+self.extensionReads.roundInput.elongation.roundDir+'/temp1'
				os.system(commendline)
		else:
			self.extensionSeqNote='onlyOneRead'
			for gseq in SeqIO.parse(self.extensionReads.extensionReads,'fasta'):
				l='>'+gseq.id+"\n"+gseq.seq+"\n"
				file1.writelines(l)
		file1.close()

	def hifiasm(self):
		self.hifiasmDir=self.extensionReads.roundInput.elongation.roundDir+"/hifiasm"
		if not os.path.exists(self.hifiasmDir):
			os.makedirs(self.hifiasmDir)
		hifiasmID=self.hifiasmDir+"/PotentialExtensionContig.asm"
		self.hifiasmResult=hifiasmID+".p_ctg.fa"
		self.hifiasmCommand="hifiasm -i -n 1 -o "+hifiasmID+" -t "+self.extensionReads.roundInput.elongation.base.thread+" -pacbio-hifi "+self.extensionReads.extensionReads
		if os.path.exists(self.hifiasmResult)!=True or os.path.getsize(self.hifiasmResult)==0:
			os.system(self.hifiasmCommand)
		
		for files in os.listdir(self.hifiasmDir):
			if '.p_ctg.gfa' in files:
				if 'bp.p_ctg.gfa' in files:
					self.hifiasmResultInput=self.hifiasmDir+"/"+files
				else:
					if 'hap' not in files:
						self.hifiasmResultInput=self.hifiasmDir+"/"+files
		
		file1=open(self.hifiasmResultInput,'r')
		file2=open(self.hifiasmResult,'w')
		hifiasmResultDict={}
		for rseq in file1:
			rseq1=rseq.rstrip().split('\t')
			if rseq1[0]=='S':
				l='>'+rseq1[1]+"\t"+'\t'.join(rseq1[3:])+"\n"+rseq1[2]+"\n"
				file2.writelines(l)
			else:
				if rseq1[1] in hifiasmResultDict:
					listhifi=hifiasmResultDict[rseq1[1]]
				else:
					listhifi=[]
				listhifi.append(rseq1[4])
				hifiasmResultDict[rseq1[1]]=listhifi
		file1.close()
		file2.close()
		self.hifiasmResultDict=hifiasmResultDict

	def readlog(self):
		logfilet=open(self.log,'r')
		self.selectContigNote=''
		for row in logfilet:
			row1=row.rstrip().split('\t')
			if row1[0]=='hifiasmDir':
				self.hifiasmDir=row1[1]
			elif row1[0]=='hifiasmCommand':
				self.hifiasmCommand=row1[1]
			elif row1[0]=='hifiasmResult':
				self.hifiasmResult=row1[1]
			elif row1[0]=='hifiasmResultDict':
				if len(row1)>1:
					row2=row1[1].split(':')
					hifiasmResultDict={}
					for i in row2:
						row3=i.split(';')
						hifiasmResultDict[row3[0]]=row3[1:]
					self.hifiasmResultDict=hifiasmResultDict
				else:	
					self.hifiasmResultDict={}

			elif row1[0]=='readCommonResult':
				self.readCommonResult=row1[1]
			elif row1[0]=='extensionSeqNote':
				self.extensionSeqNote=row1[1]
			elif row1[0]=='extensionContigs':
				self.extensionContigs=row1[1]

			elif row1[0]=='potentialExtensionContigsAln':
				self.potentialExtensionContigsAln=row1[1]
			elif row1[0]=='selectContigIdentity':
				self.selectContigIdentity=float(row1[1])
			elif row1[0]=='selectContigDistance':
				self.selectContigDistance=float(row1[1])
			elif row1[0]=='selectContigAlnLength':
				self.selectContigAlnLength=float(row1[1])
			elif row1[0]=='contigAlnMerge':
				self.contigAlnMerge=bool(row1[1])
			elif row1[0]=='contigAlnMergeIdentity':
				self.contigAlnMergeIdentity=float(row1[1])

			elif row1[0]=='potentialExtensionContigsAlnDict':
				self.potentialExtensionContigsAlnDict={}
				for i in row1[1:]:
					i1=i.split('-')
					if "ovl" in i1:
						k1,k2,k3,a,b,c,d,len1=i1
						k=k1+"-"+k2+"-"+k3
					else:
						k,a,b,c,d,len1=i1
					alist=[]
					a1=a.rstrip().split(';')
					for a2 in a1:
						if a2!='':
							a3=a2.split(',')
							alist.append([int(a3[0]),int(a3[1])])
					blist=[]
					b1=b.rstrip().split(';')
					for b2 in b1:
						if b2!='':
							b3=b2.split(',')
							blist.append([int(b3[0]),int(b3[1])])
					clist=[]
					c1=c.rstrip().split(';')
					for c2 in c1:
						if c2!='':
							c3=c2.split(',')
							clist.append([int(c3[0]),int(c3[1])])
					dlist=[]
					d1=d.rstrip().split(';')
					for d2 in d1:
						if d2!='':
							d3=d2.split(',')
							dlist.append([int(d3[0]),int(d3[1])])
					self.potentialExtensionContigsAlnDict[k]=[alist,blist,clist,dlist,int(len1)]
			elif row1[0]=='selectExtensionContigsAln':
				self.selectExtensionContigsAln=[]
				line1='\t'.join(row1[1:])
				line2=line1.split(';')
				for i in line2:
					if i!='':
						i1=i.split('\t')
						b=i1[-1].split('-')[-1]
						a='\t'.join(i1[:-1])+"\t"+i1[-1].split('-')[0]
						self.selectExtensionContigsAln.append([a+"\n",int(b)])
			elif row1[0]=='maxextensionLen':
				self.maxextensionLen=int(row1[1])
			elif row1[0]=='selectContigNote':
				self.selectContigNote=self.selectContigNote+'\t'.join(row1[1:])+"\n"
			elif row1[0]=='FuzzyAln':
				self.FuzzyAln=[]
				if len(row1)>1:
					line1=row1[1].rstrip().split(';')
					for i in line1:
						i1=i.split('-')
						self.FuzzyAln.append(i1)
			elif row1[0]=='extensionContigID':
				self.extensionContigID=row1[1].split(';')
			elif row1[0]=='extensionLength':
				self.extensionLength=int(row1[1])
			elif row1[0]=='extensionSequence':
				self.extensionSequence=row1[1]
		logfilet.close()
