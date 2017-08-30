__author__="Yesesri Cherukuri"
__date__ ="Created:  July 01,1017"
#==============================================================
#### Script to extract Extract Primary Transcripts ###########
#==============================================================
from hagsc_lib import iterFASTA, writeFASTA
from hagsc_lib import SeqRecord 
from os.path import join
#==============================================================
def real_main() : 
	#chrList = ['01','02','03','04','05','06','07','08','09']
	chrList = ['10']
	for chr in chrList : 
		outfile = open("sorghum_bicolor.CodingSequence.Sb%s.fasta"%chr,'w')
		for record in iterFASTA(open("sorghum_bicolor.CodingSequence.IDrenamed.fasta",'r')) :  
			ChrId = record.id.split("|")[3].split("_")[0]
			if not (ChrId == "Chr%s"%chr) : continue
			try : 
				count = record.id.split("|")[3].split("_")[1]
				if not (count == 1) : continue
				else : 
					r = SeqRecord( id= record.id , seq= record.seq, description='' )
					writeFASTA([r], outfile)			
			except IndexError : 
			#x = record.id.split("|")[0].split(".")[2]
			#if not (x == "1") : continue
				print ChrId 
				r = SeqRecord( id= record.id , seq= record.seq, description='' )
				writeFASTA([r], outfile)
		####
	####
####
#==============================================================
if ( __name__ == '__main__' ):
    real_main()
