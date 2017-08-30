__author__="Yesesri Cherukuri"
__date__ ="Created:  July 07,2017"

#=======================================================#
##Script to Trim the overlap between adjacent clones based on
###the tiling path constructed based on Sorghum gene locations
#=======================================================#
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
def sharedExonsGroup(sharedExonList,ExonDict) : 
	startList = []
	endList  = []
	for exon in sharedExonList :
		startList.append(ExonDict[exon][0]) 
		endList.append( ExonDict[exon][1])
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
	for exonTuple in exonDict.itervalues()   : endList.append(exonTuple[1])
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
		if(False):
			RefFile          = "sugarCaneBACs.%s.fasta"%chr
			overlapFile      = open ( "overlappingPairs.%s.dat"%chr,'w')
			TrimFile         = open ("overlappingPairs.%s.TrimCordinates.dat"%chr,'w')		
			#Read in cloneOrder 
			OrderedCloneList = []
			for line in open("%s_test_joins.dat"%chr,'r') : 
				x = line.split("10000")
				for e in x : OrderedCloneList.append( e.strip() )
			####
			cloneGeneDict    = {}
			cloneLenDict     = {}
			for line in open("sortByScaff_%s_test.out"%chr,'r'): 
				if(line.count("==") > 0 ) : continue
				x = line.split(None)
				clone    = x[2]
				cloneLen = int(x[3])
				Gstart   = int(x[4])
				Gend     = int(x[5])
				gene     = x[7]
				if not clone in cloneLenDict.keys() : cloneLenDict[clone]  = cloneLen
				try             : cloneGeneDict[clone][gene]= (Gstart,Gend)
				except KeyError : 
					cloneGeneDict[clone] = {}
					cloneGeneDict[clone][gene] = (Gstart,Gend)
			####
			#stage1 , find overlaps and Trim corordinates
			trimCloneSet     =  set()
			DumpcloneSet     =  set()
			for line in open ("DumpClone.%s.dat"%chr,'r') : DumpcloneSet.add(line.strip())
			for n in range (0,len(OrderedCloneList)-1): 
				cloneA = OrderedCloneList[n].replace("rc","").strip()
				cloneB = OrderedCloneList[n+1].replace("rc","").strip()
				if(cloneA in DumpcloneSet or cloneB in DumpcloneSet) : continue
				sharedList = list( set(cloneGeneDict[cloneA].keys()) & set(cloneGeneDict[cloneB].keys()) )
				if(len(sharedList) > 0) :
					#which clone to trim
					#add 1 for NO
					trimA = 0
					trimB = 0
					#stats on cloneA
					s1,e1    =sharedExonsGroup(sharedList,cloneGeneDict[cloneA])
					L1       = cloneLenDict[cloneA]
					TrimLen1 = min(e1,(L1-s1))
					Frag1    = L1 - TrimLen1
					Ts1,Te1  = TrimCordinates (s1,e1,L1)
					trimA    = ExoninTrimRegion(cloneGeneDict[cloneA],Ts1,Te1,trimA)
					#stats on cloneB
					s2,e2  = sharedExonsGroup(sharedList,cloneGeneDict[cloneB])
					L2     = cloneLenDict[cloneB]
					TrimLen2 = min(e2,(L2-s2))
					Frag2    = L2 - TrimLen2
					Ts2,Te2  = TrimCordinates (s2,e2,L2)
					trimB    = ExoninTrimRegion(cloneGeneDict[cloneB],Ts2,Te2,trimB)
					#Cases to decide which clone to trim
					#cloneA trim ; retention frag > cloneB trim ; Retention Frag
					if(cloneA in trimCloneSet) : 
						TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1))
						overlapFile.write("%s\t%s\n"%(cloneA,cloneB))
						continue
					elif(cloneB in trimCloneSet) : 
						TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2))	
						overlapFile.write("%s\t%s\n"%(cloneA,cloneB))			
						continue				
					elif( (Frag2+TrimLen2+Frag1) > (Frag1+TrimLen1+Frag2)  and  trimA == 0 or  trimA < trimB or  trimA == trimB) : 
						trimCloneSet.add(cloneA)
						TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneB,cloneA,Ts1,Te1))
						overlapFile.write("%s\t%s\n"%(cloneA,cloneB))
						continue
					elif ( trimB ==0 or trimB < trimA or trimA == trimB )   : 
						trimCloneSet.add(cloneB)
						TrimFile.write("%s\t%s\t%d\t%d\n"%(cloneA,cloneB,Ts2,Te2))	
						overlapFile.write("%s\t%s\n"%(cloneA,cloneB))	
						continue	
					else : assert False	
					####				
			####	
			overlapFile.close()
			TrimFile.close()
			####
			#Stage2 : Trim Clones
			outfile  = open("sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr,'w')
			mainDict = {}
			for line in open ("overlappingPairs.%s.TrimCordinates.dat"%chr,'r'): 
				x     = line.split(None)
				clone = x[1]
				start = int(x[2])
				end   = int(x[3])
				try              : mainDict[clone].append((start,end))
				except  KeyError : mainDict[clone] = [(start,end)]
			####
			cloneSeqDict = cloneSeq(RefFile)
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
			##Stage 3 :  Final Clone File
			outfile = open("sugarCaneBACs.%s.Trimmed.fasta"%chr,'w')
			TrimcloneSeqDict = cloneSeq("sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr)
			for clone in OrderedCloneList : 
				Id = clone.replace("rc","").strip()
				if(Id in DumpcloneSet): continue
				if(Id in TrimcloneSeqDict.keys()) : 
					if(clone.count("rc")>0) : 
						tmpseq = "%s"%(TrimcloneSeqDict[Id])
						newSeq = tmpseq.translate(trans_table)[::-1]
					else : newSeq = TrimcloneSeqDict[Id]
					if( len(newSeq) < 10000) : continue 
					r = SeqRecord( id= Id , seq= newSeq, description='' )
					writeFASTA([r], outfile)
					continue
				####
				else : 
					if(clone.count("rc")>0) : 
						tmpseq = "%s"%(cloneSeqDict[Id])
						newSeq = tmpseq.translate(trans_table)[::-1]	
					else : newSeq = cloneSeqDict[Id]
					r = SeqRecord( id= Id , seq= newSeq, description='' )
					writeFASTA([r], outfile)
					continue
				####
			####
			outfile.close()	
			system("rm -f sugarCaneBACs.%s.overlapTrimmed.tmp.fasta"%chr)
		##########Find the Freq of Gene ocurance in Trimmed Vs normal fasta file ###############
		if(True): 
			outfile = open("%s.geneFreq.PostDump.dat"%chr,'w')
			GeneCloneDict     = {}
			for line in open("sortByScaff_%s_test.out"%chr,'r'): 
				if(line.count("==")>0) : continue
				x = line.split(None)
				clone    = x[2]
				gene     = x[7]
				try             : GeneCloneDict[gene].append(clone)
				except KeyError : GeneCloneDict[gene] = [ clone ] 
			####
			BesthitDict     = {}
			TrimBesthitDict = {}
			#   sugarCaneBACs.sorghum_bicolor.CodingSequence.Sb01.Trimmed.New.Besthit
			for line in open("sugarCaneBACs.sorghum_bicolor.CodingSequence.%s.Besthit"%chr,'r') : 
				x = line.split(None)
				gene  = x[0]
				clone = x[5]
				ID    = float(x[10])
				cov   = float(x[11])
				try             : BesthitDict[clone][gene] = (ID,cov) 
				except KeyError : 
					BesthitDict[clone] = {}
					BesthitDict[clone][gene] = (ID,cov)
			for line in open("/home/t4c1/WORK/YESESRI_WORK/sugarcaneBACs/singleTilingPath/sugarCaneBACs.%s.PostDumpClones.Besthit"%chr,'r') : 
				x     = line.split(None)
				gene  = x[0]
				clone = x[5]
				ID    = float(x[10])
				cov   = float(x[11])
				try             : TrimBesthitDict[clone][gene] = (ID,cov) 
				except KeyError : 
					TrimBesthitDict[clone] = {}
					TrimBesthitDict[clone][gene] = (ID,cov)
			#for gene in shareGeneList : 
			TrimGeneCloneDict = {}
			for gene in GeneCloneDict.keys(): 
				for clone in GeneCloneDict[gene] : 
					try : 
						orgID   = float ( BesthitDict[clone][gene][0] )
						orgCov  = float ( BesthitDict[clone][gene][1] )
						if not (orgID>85 or  orgCov>35) : continue
						TrimId  = float ( TrimBesthitDict[clone][gene][0] )
						TrimCov = float ( TrimBesthitDict[clone][gene][1] )
					except KeyError : continue
					if not (TrimCov >= orgCov ) : continue
					print gene , clone , orgCov , TrimCov
					try             : TrimGeneCloneDict[gene].append(clone)
					except KeyError : TrimGeneCloneDict[gene] = [ clone ] 
				####
			####
			for gene in TrimGeneCloneDict.keys() : 
				if(len(TrimGeneCloneDict[gene]) > len(GeneCloneDict[gene]) ) : continue
 				outfile.write("%s\t%d\t%d\n"%(gene , len(GeneCloneDict[gene]) , len(TrimGeneCloneDict[gene])))
 			####
			outfile.close()
		####
	####Chromosome end loop
#####
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
