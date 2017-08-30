from hagsc_lib import iterFASTA, writeFASTA
from hagsc_lib import SeqRecord
from sys import argv
import os
#=======================================================================================================
def real_main():	
	basePath    = "/global/projectb/scratch/sri_ch/Para_AmbiBaseFixing"  
	#Imput Files
	##List of Ambi Basedin Ref File that can be modified
	AmbigiousBasesSummary    = "%s/AmbigiousBases.Summary"%(basePath)
	RefAmbigiousBasesSummary = "%s/AmbigiousBases.maingenome.dat"%(basePath)
	mainGenome               = "%s/Paraphysomonas_imperforata.mainGenome.scaffolds.fasta"%(basePath)
	#outfiles 
	UnmodifiedBases         = open("%s/UnModifiedBases.dat"%(basePath),'w')
	logFile	                = open("%s/Modification.log"%(basePath),'w')
	SummaryFile             = open("%s/PostModified.Summary"%(basePath),'w')
	if os.path.isfile("%s.Modified"%mainGenome): print "File *AmbiBaseModified.fasta exist"
	else : outFastaFile = open("%s.Modified"%mainGenome,'w')
	Finalcount = 0
	if(True):
		#All Ambigious Bases in Ref Genome 
		RefAmbiBaseList = []
		scaffID         = ""
		RefAmbiBaseDict = {}
		for line in open ("%s"%(RefAmbigiousBasesSummary),'r'): 
			if(line.count("scaffold")==1) : scaffID = line.replace("-","").strip()
			else : 
				Pos    = int(line.split(None)[0])
				Base   = line.split(None)[1].strip()
				makeId = "%s:%i"%(scaffID,Pos)
				RefAmbiBaseList.append(makeId)
				RefAmbiBaseDict[makeId] = Base
		####
		AmbiPosDict = {} 
		RefIdList	= []
		AmbiPosList = []
		for line in open("%s"%(AmbigiousBasesSummary),'r'):	 
			if(line.count("#") > 0 ) : continue		    
			splitline		  = line.split()
			RefId			  = splitline[0].strip()
			#0-based
			pos			      = int(splitline[1])
			RefBase		      = splitline[2]
			ConsensusReadBase = splitline[3]
			RefIdList.append(RefId)
			makeId = "%s:%i"%(RefId,pos)
			AmbiPosList.append(makeId)
			####
			try : 
				AmbiPosDict[RefId].append( (pos,RefBase,ConsensusReadBase) )
			except KeyError : 
				AmbiPosDict[RefId] = [(pos,RefBase,ConsensusReadBase)]
		####
		#Unmodified bases List
		for id in RefAmbiBaseList : 
			if not (id in AmbiPosList ) :  UnmodifiedBases.write("%s\t%i\t%s\n"%(id.split(":")[0],int(id.split(":")[1]),RefAmbiBaseDict[id]))
		####
		for	 record in iterFASTA(open("%s"%(mainGenome),'r')): 
			if not (record.id in RefIdList) : 
				 writeFASTA([record],outFastaFile)
				 continue
			####					  
			seq		  = record.seq
			logFile.write( ">%s Total bases Modified : %i\n"%(record.id,len(set(AmbiPosDict[record.id]))) )
#			Finalcount+=int(len(set(AmbiPosDict[record.id])))
			for eachTuple in sorted( set(AmbiPosDict[record.id])  ,key=lambda tup: tup[0],reverse = True): 
				#Check if Ref base in seq and reported are same
				 if(seq[int(eachTuple[0])] == eachTuple[1]) : 
					logFile.write("%s\t%s\t%s\n"%(record.id,int(eachTuple[0]),eachTuple[1].strip()))
					Finalcount+=1
					front	 = seq[:int(eachTuple[0])]
					end	 = seq[int(eachTuple[0]):]
					replace = str(eachTuple[2])
					replace+=end[1:]
					front+=replace	 
					seq	 = front			  
			####
			if(len(record.seq) == len(seq)) : 
				 new_record = SeqRecord(id = "%s"%(record.id) , seq = seq , description = '')
				 writeFASTA([new_record],outFastaFile)
			else : print "Seq Length different after modified %s"%(record.id)
		####
	SummaryFile.write("Total Ambigious Bases in Ref Genome  : %i\n"%(len(RefAmbiBaseList)))
	SummaryFile.write("Total Bases Modified                 : %i\n"%Finalcount)
	SummaryFile.write("Total Bases Not Modified             : %i\n"%(len(RefAmbiBaseList) - Finalcount))
	outFastaFile.close()
	logFile.close()
####	
#=======================================================================================================
if ( __name__ == '__main__' ):
	real_main()					   
