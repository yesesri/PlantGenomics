__author__="Yesesri Cherukuri"
__date__ ="Created:  July 07,1017"

from hagsc_lib import iterFASTA, writeFASTA
from hagsc_lib import SeqRecord 
from os.path import join
from os import system
import subprocess
from string import maketrans
trans_table = maketrans("ATGC","TACG")
#==============================================================#
def iterator(handle):
    while True:
        id   = handle.readline()
        c1   = handle.readline()
        c2   = handle.readline()
        if (not id ): return
        if (id.count(">") > 0 ): yield {'ID' : id, 'c1' : c1 , 'c2' : c2 }
        ####
    ####
#==============================================================
def sharedExonsGroup(geneOrdinateList,gene) : 
	startList = []
	endList  = []
	for n in range(0,len(geneOrdinateList)) :
		startList.append(geneOrdinateList[n][0]) 
		endList.append( geneOrdinateList[n][1])
	####
	startList.sort()
	endList.sort()
	return(startList[0],endList[-1])
#==============================================================
def TrimCordinates (s,e,L) :
	if(e < (L-s)) : return(0,e)
	else          : return(s,L)
	####
#=============================================================
def ExoninTrimRegion(exonDict,Ts,Te,trimcount) : 
	endList  = []
	for eList in exonDict.itervalues()   :
		for n in range (0,len(eList)): 
			endList.append(eList[n][1])
		####
	####
	endList.sort(reverse=True)
	for e in endList : 
		if (e>Ts and e<Te): trimcount+=1
	###
	return(trimcount)
#=============================================================
def cloneSeq(cloneFile): 
	cloneSeqDict = {}
	for record in iterFASTA(open("%s"%(cloneFile),'r')) :cloneSeqDict[record.id] = record.seq
	####
	return (cloneSeqDict)
####	
#==============================================================
def real_main() : 
	chrList          = ['Sb01','Sb02','Sb03','Sb04','Sb05','Sb06','Sb07','Sb08','Sb09','Sb10']
	for chr in chrList : 
		#if not (chr == "Sb01"): continue
		TrimFile         = open ("overlappingPairs.%s.TrimCordinates.twoCopy.dat"%chr,'w')		
		#Read in cloneOrder 
		cloneGeneDict    = {}
		cloneLenDict     = {}
		for line in open("/home/t4c1/WORK/YESESRI_WORK/sugarcaneBACs/singleTilingPath/TrimmedFiles/sugarCaneBACs.sorghum_bicolor.CodingSequence.%s.Trimmed.New.Besthit"%chr,'r'): 
			x = line.split(None)
			gene     = x[0]
			clone    = x[5]
			cloneLen = int(x[6])
			Gstart   = int(x[7])
			Gend     = int(x[8])
			if not clone in cloneLenDict.keys() : cloneLenDict[clone]  = cloneLen
			try  : 
				try             : cloneGeneDict[clone][gene].append( (Gstart,Gend) )
				except KeyError : cloneGeneDict[clone][gene] = [ (Gstart,Gend) ] 
			except KeyError: 
				cloneGeneDict[clone] = {}	
				cloneGeneDict[clone][gene] = [(Gstart,Gend) ]			
		####
		#stage1 , find overlaps and Trim corordinates
		selectcloneList = []
		trimCloneSet     =  set()
		DumpcloneSet     =  set()
		for line in open ("DumpClone.%s.dat"%chr,'r') : DumpcloneSet.add(line.strip())
		for record in iterator(open("%s.twocopyGenes.INFO.dat"%chr,'r')):	
			sharedgene  =  record['ID'].replace(">",'').strip()
			cloneA      =  record['c1'].split(":")[1].strip().replace("rc","") 
			cloneB      =  record['c2'].split(":")[1].strip().replace("rc","")	
			selectcloneList = []		 
			if(cloneA in DumpcloneSet or cloneB in DumpcloneSet) : continue
			#stats on cloneA
			print ">" , sharedgene
			print cloneA , cloneGeneDict[cloneA][sharedgene] , cloneLenDict[cloneA]
			print cloneB , cloneGeneDict[cloneB][sharedgene] , cloneLenDict[cloneB]
			#add 1 for NO
			trimA = 0
			trimB = 0
			#stats on cloneA
			s1,e1    = sharedExonsGroup(cloneGeneDict[cloneA][sharedgene] , sharedgene)
			L1       = cloneLenDict[cloneA]
			TrimLen1 = min(e1,(L1-s1))
			Frag1    = L1 - TrimLen1
			Ts1,Te1  = TrimCordinates (s1,e1,L1)
			trimA    = ExoninTrimRegion(cloneGeneDict[cloneA],Ts1,Te1,trimA)
			#stats on cloneB
			s2,e2  = sharedExonsGroup(cloneGeneDict[cloneB][sharedgene] , sharedgene)
			L2     = cloneLenDict[cloneB]
			TrimLen2 = min(e2,(L2-s2))
			Frag2    = L2 - TrimLen2
			Ts2,Te2  = TrimCordinates (s2,e2,L2)
			trimB    = ExoninTrimRegion(cloneGeneDict[cloneB],Ts2,Te2,trimB)
			#which clone to trim
			if(cloneA in trimCloneSet) : 
				TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1))
				print "%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1)
				continue
			elif(cloneB in trimCloneSet) : 
				TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2))		
				print 	"%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2)
				continue				
			elif( (Frag2+TrimLen2+Frag1) > (Frag1+TrimLen1+Frag2)  and  trimA == 0 or  trimA < trimB or  trimA == trimB) : 
				trimCloneSet.add(cloneA)
				TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1))
				print "%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1)
				continue
			elif ( trimB ==0 or trimB < trimA or trimA == trimB )   : 
				trimCloneSet.add(cloneB)
				TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2))	
				print 	"%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2)
				continue	
			else : assert False	
			####				
		####	
		TrimFile.close()
		####
		#Stage2 : Trim Clones
		outfile  = open("sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr,'w')
		mainDict = {}
		for line in open ("overlappingPairs.%s.TrimCordinates.twoCopy.dat"%chr,'r'): 
			x     = line.split(None)
			clone = x[1]
			start = int(x[2])
			end   = int(x[3])
			try              : mainDict[clone].append((start,end))
			except  KeyError : mainDict[clone] = [(start,end)]
		####
		cloneSeqDict = cloneSeq("sugarCaneBACs.%s.Trimmed.fasta"%chr)
		for clone,ordinatesList in mainDict.iteritems() :
			ordinatesList.sort(key = lambda x:x[0])
			if(len(ordinatesList)>1):  
				keepSeq_start = int(ordinatesList[0][1])
				keepSeq_end   = int(ordinatesList[-1][0])
				newSeq        = cloneSeqDict[clone][keepSeq_start:keepSeq_end]
			else : 
				#Front Trim
				if(int(ordinatesList[0][0]) == 0) : 
					keepSeq_start = int(ordinatesList[0][1])
					keepSeq_end   = None
					newSeq        = cloneSeqDict[clone][keepSeq_start:]
				#Rear Trim
				else : 
					keepSeq_start = None
					keepSeq_end = int(ordinatesList[0][0])
					newSeq 	    = cloneSeqDict[clone][:keepSeq_end]
			####
			r = SeqRecord( id= clone , seq= newSeq, description='' )
			writeFASTA([r], outfile)	
		####
		outfile.close()
		##Stage3 : Final Clone File
		outfile = open("sugarCaneBACs.%s.Trimmed.final.fasta"%chr,'w')
		TrimcloneSeqDict = cloneSeq("sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr)
		for Id in cloneSeqDict.keys() : 
			if(Id in DumpcloneSet): continue
			if(Id in TrimcloneSeqDict.keys()) : 
				newSeq = TrimcloneSeqDict[Id]
				if( len(newSeq) < 10000) : continue 
				r = SeqRecord( id= Id , seq= newSeq, description='' )
				writeFASTA([r], outfile)
				continue
			else : 
				newSeq = cloneSeqDict[Id]
				r = SeqRecord( id= Id , seq= newSeq, description='' )
				writeFASTA([r], outfile)
				continue
			####
		####
		outfile.close()	
		system("rm -f sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr)
	####Chromosome end loop
#####
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
