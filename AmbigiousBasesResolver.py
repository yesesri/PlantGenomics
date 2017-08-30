__author__="ycherukuri@hudsonalpha.org"
__date__ ="Created:	 02/24/2017"

from hagsc_lib import iterFASTA, writeFASTA
from sys import argv
import subprocess 
import re
from os.path import join
from sys import stderr
#===============================================================
def LocateAmbiBasesInGenome(FastaFile,RefIdList): 
	#List of non ambigious bases
	BaseList				  = ["A","T","G","C","a","t","g","c","N","n"]
	#Intilization
	##positions of amibgious bases on each scaffold	 
	RefAmbiBaseLocDict		  = {}
	##50bp Flanking REF seq on either side of the Ambigious base 
	RefSeqDict				  = {}
	TotalAmbiBaseCount		  = 0
	stderr.write( "Scanning %s for Ambigious Bases\n"%FastaFile)
	# iterate Reference Genome file
	for record in iterFASTA(open("%s"%FastaFile)) : 
		# if Ref scaffold has No Alignments ( from bam parse) do nothing and continue to next scaffold 
		if not (record.id in RefIdList) : continue		  
		#0-based position 
		#current base position
		baseIndex = 0
		for base in record.seq : 
			if(base in BaseList): baseIndex+=1
			#Gather the positions of amibgious bases on each scaffold				
			else : 
				try				  : RefAmbiBaseLocDict[record.id].append((baseIndex,base))
				except KeyError	  : RefAmbiBaseLocDict[record.id] = [(baseIndex,base)]			
				#Pull 50bp flanking ref seq on either side of the ambigious base		
				key	   = "%s:%s"%(record.id,baseIndex)
				tmpSeq =   "%s%s"%(record.seq[(baseIndex - 50 ) : baseIndex] , record.seq[ baseIndex :(baseIndex + 50 )] )
				RefSeqDict[key] = ( (tmpSeq,base ) ) 
				#next base position
				baseIndex+=1
			  ####
		 #### base loop of each read seq
	#### end loop of fastafile
	return(RefAmbiBaseLocDict,RefSeqDict)
####
#===============================================================
def ParseBamAlignments(bamfile):
	AlignSeqDict		= {}
	AlignStartPosDict	= {}  
	RefIdList			= [] 
	cigarRE				= re.compile(r'(\d+)([A-Z]{1})').findall
	#Read in alignments from bamfile 
	stderr.write("Parsing %s\n"%bamfile)
	cmd = "samtools view -sort %s"%(bamfile)
	#cmd = "samtools view -S  -F 4	%s"%bamfile 
	p	= subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout :	
		splitline			 = line.split()
		LibReadID			 = splitline[0]
		RefID				 = splitline[2]
		#convert 1-based position to 0-based position
		alignStartPos		 = int(splitline[3]) - 1
		QueryReadAlignSeq    = splitline[9]
		cigarList			 = cigarRE(splitline[5]) 
		Len = 0
		#Modify Align Reads based on CIGAR string
		#######################################
		##SoftClip		 : Trim the Sequence
		##HardClip		 : No change 
		##Match/Mismatch : Nochange 
		##Insertion		 : Delete Base
		##Deletion		 : Insert "_" 
		#######################################
		for cigar in cigarList :
			string_1 = ""
			string_2 = ""
			string_3 = ""
			if(cigar[1]=='D'): 
				string_1 = QueryReadAlignSeq[:Len]
				string_2 = int(cigar[0])*"_"
				string_3 = QueryReadAlignSeq[Len:]
				string_1+=string_2
				string_1+=string_3
				QueryReadAlignSeq = string_1
			if(cigar[1]=='S'):
				string_1 = QueryReadAlignSeq[:Len]
				string_2 = QueryReadAlignSeq[Len+int(cigar[0]):]
				string_1+=string_2
				QueryReadAlignSeq = string_1			
			if(cigar[1] == 'I'):
				string_1 = QueryReadAlignSeq[:Len]
				string_2 = QueryReadAlignSeq[Len+int(cigar[0]):]
				string_1+=string_2
				QueryReadAlignSeq = string_1					
			####			
			if not (cigar[1]=='S' or cigar[1]=='H' or cigar[1]=='I') : Len+=int(cigar[0])	
		####
		#check if everything went well
		if not (len(QueryReadAlignSeq) == Len) : print "Someting Went Wrong"		
		####
		# List out all the Read Alignment start Positions on each REF scaffold				
		try : 
			AlignStartPosDict[RefID].append( alignStartPos )
		except KeyError : 
			AlignStartPosDict[RefID] = [ alignStartPos ]
			#List of all REFERENCE Scaffolds  having read alignments
			RefIdList.append(RefID)
		####
		makeID			 = "%s:%i"%(RefID,alignStartPos) 
		# gather all the Reads aligned to REF scaffold with same alignStart
		try : 
			AlignSeqDict[makeID].append( (LibReadID , QueryReadAlignSeq) )
		except KeyError : 
			AlignSeqDict[makeID] = [(LibReadID , QueryReadAlignSeq)]			
		####			
	#### bamfile end loop
	p.poll() 
	return( AlignStartPosSorter(AlignStartPosDict),AlignSeqDict,RefIdList)	 
#### end of function			
#======================================================================
def AlignStartPosSorter(AlignStartPosDict) : 
	NewAlignStartPosDict = {}
	for RefId,AlignStartPosList	 in	 AlignStartPosDict.iteritems() : 
		 NewAlignStartPosDict[RefId]	=	 sorted(set(AlignStartPosList)) 
	return ( NewAlignStartPosDict ) 
#### end of function
#=================================================================
def ConBaseTester(RefBase,ConBase,BaseCountDict):
	if(RefBase == "R"): 
		 if( ConBase == "A" or ConBase == "G") : return(RefBase,ConBase)
		 elif(int(BaseCountDict['A']) > int(BaseCountDict['G']) ) : 
			 ConBase = "A"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "G"
			 return(RefBase,ConBase) 
	if(RefBase == "Y") : 
		 if( ConBase == "C" or ConBase == "T") : return(RefBase,ConBase)
		 elif(int(BaseCountDict['C']) > int(BaseCountDict['T']) ) : 
			 ConBase = "C"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "T"
			 return(RefBase,ConBase) 
	if(RefBase == "S") : 
		 if( ConBase == "G" or ConBase == "C")	: return(RefBase,ConBase)
		 elif(int(BaseCountDict['G']) > int(BaseCountDict['C']) ) : 
			 ConBase = "G"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "C"
			 return(RefBase,ConBase) 			   
	if(RefBase == "W" ) :
		 if(ConBase == "A" or ConBase == "T") : return(RefBase,ConBase)
		 elif(int(BaseCountDict['A']) > int(BaseCountDict['T']) ) : 
			 ConBase = "A"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "T"
			 return(RefBase,ConBase) 
	if(RefBase == "K"): 
		 if( ConBase == "G" or ConBase == "T")	: return(RefBase,ConBase)
		 elif(int(BaseCountDict['G']) > int(BaseCountDict['T']) ) : 
			 ConBase = "G"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "T"
			 return(RefBase,ConBase) 
	if(RefBase == "M" ) : 
		 if(ConBase == "A" or ConBase == "C") : return(RefBase,ConBase)
		 elif(int(BaseCountDict['A']) > int(BaseCountDict['C']) ) : 
			 ConBase = "A"
			 return(RefBase,ConBase)
		 else :	  
			 ConBase = "C"
			 return(RefBase,ConBase) 			
	if(RefBase == "B"): 
		 tmpList = []
		 if( ConBase == "C" or ConBase == "G" or ConBase == 'T') : return(RefBase,ConBase)
		 else:
			 tmpList.append(int(BaseCountDict['C']))
			 tmpList.append(int(BaseCountDict['G']))
			 tmpList.append(int(BaseCountDict['T']))
			 sorted(tmpList,reverse = True)
			 ConBase = tmpList[0]
			 return(RefBase,ConBase)
	if(RefBase == "D") :
		 if(ConBase == "A" or ConBase == "G" or ConBase == 'T') : return(RefBase,ConBase)
		 else:
			 tmpList.append(int(BaseCountDict['A']))
			 tmpList.append(int(BaseCountDict['G']))
			 tmpList.append(int(BaseCountDict['T']))
			 sorted(tmpList,reverse = True)
			 ConBase = tmpList[0]
			 return(RefBase,ConBase)
	if(RefBase == "H") :
		 if(ConBase == "C" or ConBase == "A" or ConBase == 'T') : return(RefBase,ConBase)
		 else:
			 tmpList.append(int(BaseCountDict['C']))
			 tmpList.append(int(BaseCountDict['A']))
			 tmpList.append(int(BaseCountDict['T']))
			 sorted(tmpList,reverse = True)
			 ConBase = tmpList[0]
			 return(RefBase,ConBase)			
	if(RefBase == "V" ) :
		 if(ConBase == "C" or ConBase == "G" or ConBase == 'A') : return(RefBase,ConBase)
		 else:
			 tmpList.append(int(BaseCountDict['A']))
			 tmpList.append(int(BaseCountDict['C']))
			 tmpList.append(int(BaseCountDict['G']))
			 sorted(tmpList,reverse = True)
			 ConBase = tmpList[0]
			 return(RefBase,ConBase)					   
#=================================================================
def real_main():
	basePath		  =	 "/global/projectb/scratch/sri_ch/Para_AmbiBaseFixing"
	#InputFiles 
	mainGenomeFile	  =	 "%s/Paraphysomonas_imperforata.mainGenome.scaffolds.fasta"%(basePath)	  
	bamfile			  =	 "%s/Para_SNPcall/Para.GRP_%s.gatk.bam"%(basePath,argv[1])
	##outfiles
	summaryFile		  =	 open("%s/outfiles/AmbigiousBases.%s.Summary"%(basePath,argv[1]),'w')
	GenomeAmbiInfo	  =	 open("%s/outfiles/AmbigiousBases.maingenome.%s.dat"%(basePath,argv[1]), 'w')
	MainFile		  =	 open("%s/outfiles/AmbigiousBases.%s.main"%(basePath,argv[1]), 'w')	
	#####################################################################################	
	##### DONT FORGET TO CHANGE OUTPUT FORMATTING BASED ON LEN OF READ ID ###############  
	#####################################################################################
	#step 1 : Extract Reads Aligning to each Ref scaffold
	AlignStartPosDict,AlignSeqDict,RefIdList	= ParseBamAlignments(bamfile)
	#step 2 : Locate the Ambigious Bases
	RefAmbiBaseLocDict,RefSeqDict				= LocateAmbiBasesInGenome("%s"%join(basePath,mainGenomeFile),RefIdList)
	#Step 3 : iterate Dictionaries from step 1 & 2 and Map Read to Ref Scaffold
	TotalAmbiBaseCount			= 0
	AmbiBaseAlignSeqDict		= {}
	#iterate the List of Ambigious bases on each REF scaffold
	for ReadId,AmbiList in RefAmbiBaseLocDict.iteritems() : 
		GenomeAmbiInfo.write( "-------%s--------\n"%ReadId)
		for index,base in AmbiList : 
			GenomeAmbiInfo.write("%i\t%s\n"%(index,base))
			TotalAmbiBaseCount+=1
			indexalignSeqList = []
			key = "%s:%s"%(ReadId,index)
			for AlignStartPos in AlignStartPosDict[ReadId] : 
				k = "%s:%s"%(ReadId,AlignStartPos)
				for seqID , Seq in AlignSeqDict[k] : 
					AlignEndPos	 =	int( AlignStartPos) + len(Seq)
					if(int(index)> AlignStartPos and int(index) < AlignEndPos ):
						lowerbound	 =	int(index)		 - int(AlignStartPos)
						upperbound	 =	int(AlignEndPos) -	int(index)	 
						left  = ""
						right = ""			  
						if(lowerbound > 50 ): left	  = Seq[(lowerbound - 50 ) : lowerbound]
						elif(lowerbound ==0): 
							Tmpleft	   = Seq[lowerbound]
							left  = "-"*49							  
							left+=Tmpleft																												   
						else: 
							Tmpleft	   = Seq[:lowerbound]
							left	   = "_"*(50 - len(Tmpleft))
							left+=Tmpleft
						####						
						if(upperbound > 50 )		 : right   = Seq[lowerbound: (lowerbound+ 50) ]
						elif(upperbound == len(Seq)) : 
							right		  = Seq[lowerbound] 
							Tmpright	  = "_"*49
							right+=Tmpright																					 
						else:
							right		= Seq[lowerbound:]							 
							Tmpright	= "_"*(50 - len(right))
							right+=Tmpright									 
						#print left ,":", right																			 
						left+=right
						tmpbase = Seq[lowerbound]
						indexalignSeqList.append((seqID,left,tmpbase))
					else : continue
				#### Seq at each AlignStart Pos, end loop 
			#### Reads Aligned to Ref Read, end loop  
			AmbiBaseAlignSeqDict[key] = indexalignSeqList 
		#### Ref Read Ambi Positions, end loop	
	#### Ref Read, end loop				  
	GenomeAmbiInfo.close()
	stderr.write("#AmbigiousBases in Group%s %s\n"%(argv[1],TotalAmbiBaseCount))
	MainFile.write( "Total Number of Amibiguous Bases in  %s groupfile %s file : %s \n "%(mainGenomeFile,argv[1],TotalAmbiBaseCount))
	#Step 4: Consensus Base Calling 
	stderr.write("ConsensusCalling\n")
	for key in RefSeqDict.iterkeys(): 
		 MainFile.write("====================================================================================================\n")	
		 MainFile.write("> Read ID : %s \t RefPosition(0-based) : %s \t Base  %s\n"%(key.split(":")[0] , key.split(":")[1],RefSeqDict[key][1] ) ) 
		 #Formating Output with equal spacing
		 MainFile.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>  REFERENCE\t %s \t %s \n"%(RefSeqDict[key][1],RefSeqDict[key][0]))
		 #Consensus Base calling 
		 BaseCountDict = {'A' : 0 , 'T' : 0 , 'G' :0 , 'C' :0 ,'N' :0 , 'a':0 ,'t':0 ,'g' : 0 , 'c' : 0 , 'n' :0,'_' : 0} 
		 for AlignSeq in AmbiBaseAlignSeqDict[key]: 
			 BaseCountDict[AlignSeq[2]]+=1
			 #Formating Output with equal spacing 
			 if(len(AlignSeq[0])<40): MainFile.write("%s  \t %s \t %s\n"%(AlignSeq[0],AlignSeq[2],AlignSeq[1]))
			 else: MainFile.write("%s \t %s \t %s\n"%(AlignSeq[0],AlignSeq[2],AlignSeq[1]))	 
		 tmpList = []	   
		# del BaseCountDict['_']		 
		# for values  in BaseCountDict.itervalues() : tmpList.append(values)
		# tmpList.sort(reverse = True)
		 maxcount = 0
		 del BaseCountDict['_']
		 for k,v in sorted(BaseCountDict.items(), key=lambda s: s[0]) : 
			 if(v>maxcount ): 
				maxcount	   = v 
				ConsensusBase  =  "%s"%k  
		####
		 RefBase,ConBase = ConBaseTester(RefSeqDict[key][1],ConsensusBase,BaseCountDict)	  
		#	  if(v == int(tmpList[0])) : ConsensusBase	=  "%s"%k	
		 ####
		 MainFile.write("------------------------------------------------------------------------------\n")
		 MainFile.write(" SUMMARY : \n")
		 MainFile.write("#ReadId ; Pos ; RefBase ; ConsesusBase ; depth(reads aligned at that AmbibasePos) ;#reads supporting ConsensuBase\n")
		 MainFile.write("%s\t%s\t%s\t%s\t%i\t%i\n"%(key.split(":")[0] , key.split(":")[1] ,RefBase,ConBase,len( AmbiBaseAlignSeqDict[key] ) , maxcount	))
		 MainFile.write("BaseCounts :	A  %s \t T %s \t  G %s \t C %s \t N %s \n"%(BaseCountDict['A'],BaseCountDict['T'],BaseCountDict['G'],BaseCountDict['C'],int(BaseCountDict['N'])))
		 MainFile.write("BaseCounts :	a  %s \t t %s \t  g %s \t c %s \t n %s \n"%(BaseCountDict['a'],BaseCountDict['t'],BaseCountDict['g'],BaseCountDict['c'],int(BaseCountDict['N'])))		  
		# MainFile.write("=====================================================================================================\n")
		 if(len( AmbiBaseAlignSeqDict[key] ) == 0) : continue
		 #ReadId ; Pos(0-based) ; RefBase ; ConsesusBase ; depth = #reads aligned at that AmbibasePos ; #reads supporting ConsensuBase
		 summaryFile.write("%s\t%s\t%s\t%s\t%i\t%i\n"%(key.split(":")[0] , key.split(":")[1] , RefBase , ConBase , len( AmbiBaseAlignSeqDict[key] ) , maxcount	))
	####
	stderr.write("Done!!!\n")
	MainFile.close()
	summaryFile.close()
####
#==============================================================
if ( __name__ == '__main__' ):
	real_main()
