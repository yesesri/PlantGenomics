__author__="Yesesri Cherukuri"
__date__ ="Created:  July 20 ,1017"

from hagsc_lib import iterFASTA, writeFASTA
from hagsc_lib import SeqRecord 
from os.path import join
from os import system
import subprocess
from sys import argv
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
def cloneSeq(RefFile): 
	cloneLenDict = {}
	for record in iterFASTA(open("%s"%(RefFile),'r')) :
		cloneLenDict[record.id] = len(record.seq)
	####
	return (cloneLenDict)
####	
#===============================================================
def real_main() :
	chrList  = ['Sb01','Sb02','Sb03','Sb04','Sb05','Sb06','Sb07','Sb08','Sb09','Sb10']
	for chr in chrList : 
		# > Two copy Redundant Genes , dump clones
		if(False) : 
			RefFile      	   = "sugarCaneBACs.%s.fasta"%(chr)
			sortByScaffFile    = "sortByScaff_%s_test.out"%(chr)
			geneFreqDict       = {}
			GeneClonecountDict = {}
			cloneLenDict       =  cloneSeq(RefFile)
			cloneGenecountDict = {}
			outfile            = open("DumpClone.INFO.%s.dat"%chr,'w')
			DumpFile           = open("DumpClone.%s.dat"%chr,'w')
			#Frequency of Genes (excluded genes occur multiple times on same clone)  in BAC set 
			for line in open("%s"%(sortByScaffFile),'r'):
				if(line.count("=")>0) : continue
				gene  = line.split(None)[7].strip()
				clone = line.split(None)[2].strip()
				try             : 
					geneFreqDict[gene]+=1
					GeneClonecountDict[gene].append(clone)
				except KeyError : 
					geneFreqDict[gene] = 1
					GeneClonecountDict[gene] = [ clone ] 
				####
				try             : cloneGenecountDict[clone].append(gene)
				except KeyError : cloneGenecountDict[clone] = [gene]
			####
			#Remove multiple copy genes
			DumpCloneSet = set()
			for gene,Freq in geneFreqDict.iteritems(): 
				if(Freq <3) : continue
				outfile.write(">%s\t%s\n"%(gene , Freq))
				cloneList      = GeneClonecountDict[gene]
				selecCloneList = []
				for clone in cloneList : 
					if(clone in DumpCloneSet) : continue
					selecCloneList.append( ( clone , len(cloneGenecountDict[clone]) , cloneLenDict[clone] ) )
				####
				selecCloneList.sort(key=lambda x: (x[1], x[2] ) )
				if(len(selecCloneList) < 2) : continue
				outfile.write("%s\t%s\t%d\t%d\n"%("selected",selecCloneList[-1][0],selecCloneList[-1][1],selecCloneList[-1][2]))
				for n in range(0,len(selecCloneList)-1) :
					outfile.write("%s\t%s\t%d\t%d\n"%("Dumped",selecCloneList[n][0],selecCloneList[n][1],selecCloneList[n][2]))
					DumpCloneSet.add(selecCloneList[n][0])
				outfile.write( "----------------------------------\n")
			####
			for clone in DumpCloneSet : DumpFile.write("%s\n"%clone.strip())
			outfile.close()
		####
		#Two copy Genes on non adjacent clones dump 
		if(False): 
			#if not (chr == "Sb01") : continue
			RefFile      	   = "sugarCaneBACs.%s.Trimmed.fasta"%(chr)
			sortByScaffFile    = "sortByScaff_%s_test.out"%(chr)
			GeneClonecountDict = {}
			cloneLenDict       =  cloneSeq(RefFile)
			cloneGenecountDict = {}
			outfile            = open("DumpClone.twocopyGenes.INFO.%s.dat"%chr,'w')
			DumpFile           = open("DumpClone.twocopyGenes.%s.dat"%chr,'w')
			#Frequency of Genes (excluded genes occur multiple times on same clone)  in BAC set 
			for line in open("%s"%(sortByScaffFile),'r'):
				if(line.count("=")>0) : continue
				gene  = line.split(None)[7].strip()
				clone = line.split(None)[2].strip()
				try             : GeneClonecountDict[gene].append(clone)
				except KeyError : GeneClonecountDict[gene] = [ clone ] 
				####
				try             : cloneGenecountDict[clone].append(gene)
				except KeyError : cloneGenecountDict[clone] = [gene]
			####
			DumpCloneSet   = set()
			for record in iterator(open("%s.twocopyGenes.INFO.dat"%chr,'r')):	
				gene  =  record['ID'].replace(">",'').strip()
				c1    =  record['c1'].split(":")[1].strip() 
				c2    =  record['c2'].split(":")[1].strip()
			####
				outfile.write(">%s ; %s ; %s \n"%(gene,c1,c2))
				#Remove multiple copy genes
				try : cloneList      = GeneClonecountDict[gene]
				except KeyError : 
					print gene , chr 
					continue
				selecCloneList = []
				for clone in cloneList : 
					if(clone in DumpCloneSet) : continue
					selecCloneList.append( ( clone , len(cloneGenecountDict[clone]) , cloneLenDict[clone] ) )
				####
				selecCloneList.sort(key=lambda x: (x[1], x[2] ) )
				if(len(selecCloneList) < 2) : continue
				outfile.write("%s\t%s\t%d\t%d\n"%("selected",selecCloneList[-1][0],selecCloneList[-1][1],selecCloneList[-1][2]))
				outfile.write("%s\t%s\t%d\t%d\n"%("Dumped",selecCloneList[0][0],selecCloneList[0][1],selecCloneList[0][2]))
				DumpCloneSet.add(selecCloneList[0][0])
				outfile.write( "----------------------------------\n")
			####
			for clone in DumpCloneSet : DumpFile.write("%s\n"%clone.strip())
			outfile.close()
			
			NewRefFile = open("sugarCaneBACs.%s.Trimmed.final.fasta"%chr,'w')
			for record in iterFASTA( open("%s"%RefFile ,'r') ): 
				if(record.id in DumpCloneSet) : continue
				else : writeFASTA([record], NewRefFile)
			####
			NewRefFile.close()
		####
		if(True): 
			DumpCloneSet = set()
			NewRefFile = open("sugarCaneBACs.%s.PostDumpClones.fasta"%chr,'w')
			for line in open("DumpClone.%s.dat"%chr,'r'): DumpCloneSet.add(line.strip())	
			for record in iterFASTA(open("sugarCaneBACs.%s.fasta"%chr,'r')) : 
				if record.id in DumpCloneSet : continue
				writeFASTA([record], NewRefFile)	
	####Chromosome End Loop 
####
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
