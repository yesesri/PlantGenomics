# PlantGenomics
Scripts uses internal modules for parsing FASTA/Q files, which can be replaced by Biopython modules.
Ambigious base Fixing in Post polished genomes [
      AmbigiousBasesResolver.py
        - Require bam alignment file, libraries aligned to draft genome in which ambgious bases are to be fixed 
      ModifyAmbiBases.py
        - base modified based on base depth at the position and IUPAC Neuclotides table  ].

Extract Primary Transcripts [ PrimaryTranscriptExtractor.py ]
Script to extract Haplotype Markes [ SNPmerGenerator.py ] 
Single Tiling Path Generation Scripts [ RemoveRedundantClones.py ,TrimAdjacentOverlaps.py ,TrimNonAdjOverlaps.py ]
