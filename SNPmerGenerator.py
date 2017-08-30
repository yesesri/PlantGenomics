__author__="ycherukuri@hudsonalpha.org"
__date__ ="Created:	 05/17/2017"

from hagsc_lib import iterFASTA, writeFASTA
from sys import argv
import subprocess 
import re
from os.path import join
from sys import stderr
#===============================================================
def RefLenInfo(ReferenceFile,BACcontigNames): 
	RefLenDict		= {}
	RefSeqDict      = {}
	contigNameList  = [] 
	for contigName	 in open("%s"%BACcontigNames,'r') : contigNameList.append(contigName.strip())
	for record	 in iterFASTA(open("%s"%ReferenceFile,'r')):
		if not (record.id in contigNameList ) : continue
		RefLenDict[record.id] = int(len(record.seq))
		RefSeqDict[record.id] = "%s"%record.seq		
	####
	return (RefLenDict,RefSeqDict)
####	
#===============================================================		
def	 snpInfo(snpInfofile):
	snpInfoDict = {}
	for line in open("%s"%(snpInfofile),'r'):
		lineSplit	= line.split()
		contigName = lineSplit[0].strip()
		snpPos		= int(lineSplit[1].strip()) 
		try : 
			snpInfoDict[contigName].append(snpPos)
		except KeyError : 
			snpInfoDict[contigName]= [snpPos]
		 ####
	####
	return (snpInfoDict)
####
#==============================================================
def ParseBamAlignments(LibList,contigName,Wstart,Wstop,mersize,Refmersize,RefSeqDict):	
	#Initialization
	RegionSeqDict	= {}
	merSearchDict   = {}	
	LibReadDict     = {} 
	RefmerDict      = {}
	NewRefmerDict   = {}
	makeID	        = "%s:%d-%d"%(contigName.strip(),int(Wstart)-1,int(Wstop)-1) 
#	print makeID
#	stderr.write ("Region %s\n"%(makeID))
	#Pull Reference seq within window cordinates  ; convert cordinates to 0-based
	WRefSeq         = RefSeqDict[contigName.strip()][int(Wstart)-1:int(Wstop)]
	#Extract 8 mers from the Ref Seq	
	for n in range (0,mersize,Refmersize): 
		try             : RefmerDict[WRefSeq[n:n+Refmersize]].append((n,n+Refmersize))
		except KeyError : RefmerDict[WRefSeq[n:n+Refmersize]]= [(n,n+Refmersize)]
	####
	#remove 8 mer's that occur more than once 
	for mer,FreqList in RefmerDict.iteritems(): 
		if (len(FreqList) == 1) : 
			NewRefmerDict[mer] = FreqList
			merSearchDict[mer] = {}
		else : continue
	####
	#No unique mers in the current window ; No SNP mers will be extracted 
	if(len(NewRefmerDict.keys()) == 0) : 
		stderr.write ("No Unique mer in the 80mer Ref Seq. \n No SNPmers extracted from the region : \t  %s"%makeID)
		return(None)
	else :	
		for eachLib in LibList : 	
			stderr.write("Processing %s\n"%eachLib)
			totalAlignments = 0
			#Pull the Alignments within the window cordinates 
			cmd             = 'samtools view  -q 50  /global/projectb/scratch/sri_ch/sugarcane/STP_kmerCalls/%s/%s.gatk.bam  "%s:%d-%d" '%(eachLib.strip(),eachLib.strip(),contigName,int(Wstart),int(Wstop))
#			print cmd
			p	            = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
			for line in p.stdout :	
				splitline			   = line.split()
				LibReadID			   = "%s:%i"%(splitline[0].strip(),int(splitline[1]))
				RefID				   = splitline[2]
				alignStartPos		   = int(splitline[3]) 
				totalAlignments+=1
				LibReadDict[LibReadID] = splitline[9].strip()	
				#search all the 8 mers of Ref seq in current align read sequence
				for mer in NewRefmerDict.keys() : 
				 	#Regular expression of each mer
					merRE   = re.compile(r'%s'%mer).finditer
					for eachmatch in merRE(splitline[9].strip()):
						try             : merSearchDict[mer][LibReadID].append(eachmatch.span())
						except KeyError : merSearchDict[mer][LibReadID] = [eachmatch.span()]
					#### 
				####
			####
			p.poll()
		#### end of library loop (8 libraries
	####end of  if else loop
	#mer found in all or max no of alignments
	prevAligncount = 0
	for mer,meralignDict in merSearchDict.iteritems(): 
		currentAligncount = len(meralignDict.keys())
		if(currentAligncount > prevAligncount) : 
			Finalmer       = "%s"%(mer)
			prevAligncount = currentAligncount
		####
	####
	#No Alignments within the window or Ref mer's are not found in any of the align reads
	if not  (prevAligncount == 0) : 	 	
		for LibreadID,mermatchList in merSearchDict[Finalmer].iteritems():
			Refstart = int(RefmerDict[Finalmer][0][0])
			#Ref mer occurs more than once on the read
			if(len(mermatchList)>1) : continue
			ReadmatchStart = int(mermatchList[0][0])
			if(Refstart > ReadmatchStart ) : continue
			readstart      = ReadmatchStart - Refstart
			tmpKmer        =  LibReadDict[LibreadID.strip()][readstart:]
			if(len(tmpKmer)< mersize)     : continue
			elif(len(tmpKmer) > mersize ) : SNPmer   = tmpKmer[:mersize]
			elif(len(tmpKmer)==mersize)   : SNPmer   = tmpKmer	
			try             : RegionSeqDict[SNPmer]+=1		
			except KeyError : RegionSeqDict[SNPmer] = 1		
		####
		return(RegionSeqDict,makeID,WRefSeq,Finalmer)
	else : return(None)	 		   
#### end of function 
#===============================================================			 
def real_main() : 
	mainDIct   = {}
	#Inputs [Manual Entry ]
	Refmersize           = 8
	mersize	             = 80
	LibList              = ['IGIS','IGIN_L6','IGIN_L5','IGIP_L4','IGIR','IGIP_L3','IGIQ_L2','IGIQ_L1']	
##	LibList              = ['ICGS']
	#Input File Names 
	basepath		     = "/global/projectb/scratch/sri_ch/sugarcane/STP_kmerCalls" 
	#Splitted the contigs into set of 10, so multiple jobs can be submitted
	BACcontigNames	     = "%s/splitfiles/InitialSplit/%s"%(basepath,argv[1])
##	BACcontigNames	     = "%s/Final_BAC_SET.dat"%(basepath)
	ReferenceFile	     = "%s/final.Trimmed.joined.fasta"%(basepath)
	#Outfiles
	summaryfile		     = open("%s/kmeroutfiles/%s.%s.summary.MQ50.dat"%(basepath,argv[1],argv[2]),'w')
#	summaryfile.write("#contig:region    totaldepth      #of snpmers     count of each snpmer\n")	
	InfoFile             =  open("%s/kmeroutfiles/%s.%s.INFO.MQ50.dat"%(basepath,argv[1],argv[2]),'w') 
    #pull seq and corresponding length of each Ref Contig
	RefLenDict,RefSeqDict = RefLenInfo(ReferenceFile,BACcontigNames)
	for contigName in RefLenDict.iterkeys():	
		#if not (contigName == "contig_3315"): continue				 
		# calculate upper range for each contig in Ref file 
		upperrange =  int(RefLenDict[contigName.strip()]) - int(mersize) 
		#Pull Alignments for each Window on current contig
#		for Wstart in range(1, upperrange, mersize ) : 
#			Wstop	=  Wstart + (mersize -1)
		for line in open("/global/projectb/scratch/sri_ch/sugarcane/STP_kmerCalls/kmeroutfiles/Chromosomes_failedJobs/%s"%argv[2],'r'): 
			x = line.split(":")
			Wstart = int( x[1].split("-")[0] )
			Wstop  = int( x[1].split("-")[1] )
		#	if not (Wstart ==36641) : continue
			try:
			#if(True) : 	
				RegionSeqDict,makeID,RefSeq,Finalmer  = ParseBamAlignments(LibList,contigName,Wstart,Wstop,mersize,Refmersize,RefSeqDict)  
				mainDIct[makeID]	                  = RegionSeqDict	  					
				tmpList    = []
				depth      = 0
				totaldepth = 0	
				InfoFile.write("%s\t%s\tREF\t%s\n"%(makeID,RefSeq,Finalmer))
			#	print ("%s\t%s\tREF\t%s\n"%(makeID,RefSeq,Finalmer))		
				for k,v in sorted(RegionSeqDict.iteritems(), key=lambda x: x[1] , reverse = True): 
					totaldepth+=v
					depth+=v
					tmpList.append(v)
					InfoFile.write("%s\t%s\t%i\n"%(makeID,k,v)) 
				    #if(v>1):
				    ####
			#		print("%s\t%s\t%i\n"%(makeID,k,v))
				InfoFile.write("=====================================\n")
			#	print "==============================================="
				summaryfile.write("%s\t%d\t%d\t%s\n"%(makeID,totaldepth,len(tmpList),sorted(tmpList,reverse = True))) 				
			#if (ParseBamAlignments(LibList,contigName,Wstart,Wstop,mersize,Refmersize,RefSeqDict) == None) : 
			except TypeError : 
				makeID	= "%s:%d-%d"%(contigName.strip(),int(Wstart)-1,int(Wstop)-1)
				InfoFile.write("%s\tNO ALIGNMENTS\n"%(makeID))
			#	print("%s\tNO ALIGNMENTS\n"%(makeID))				
				InfoFile.write("=====================================\n")
			#	print "==============================================="						
				summaryfile.write("%s\t0\t0\t[]\n"%(makeID))			
	    #### window end loop 
	#### Contig end loop
	InfoFile.close()
	summaryfile.close()
####
#==============================================================
if ( __name__ == '__main__' ):
	real_main()
