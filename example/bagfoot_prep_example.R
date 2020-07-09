# example
# This script shows how to prepare required data files to run "BigFoot"
# 
# By Songjoon Baek

library('bagfoot');

bamfile1= 'example_data/GH1047_PE_sorted.bam';  

cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.
countCutcountOnSites <- function(cutcount, sites)
cutcountfile1 = makeCutCountBAM(bamfile1);   # generate a BedGraph file with DNase Cleavages counts

MAPPABILITY_FILES_DIRECTORY_HG19<<-'~/Data/hg19/35mers/';   # directory of mappability files for the human hg19 ref. genome ( chr1b.out, chr2b.out, etc..)
MAPPABILITY_FILES_DIRECTORY_MM9<<-'~/Data/mm9/35mers/';     # directory of mappability files for the mouse mm9 ref. genome ( chr1b.out, chr2b.out, etc..)
 
hotspotfile1 = 'example_data/fed1_10000_hotspot.csv';   # A peak sites file is comma-separated and must has headers with chromosome, start, end (1-based)
hotspotfile2 = 'example_data/fasted1_10000_hotspot.csv';

combinedhotspotfile = combineTwoHotspots( hotspotfile1, hotspotfile2, 'fed1', 'fasted1');  # Hotspots in hotspotfile1 and hotspotfile2 are combined and saved to a new file name

tabNoMappability = MakeBiasCorrectionTableBAM(        # This function call generates a hexamer bias frequency table from the given bamfile without mappability assumption
	bamfile=bamfile1,
	outfile="Hexamer_fed_mm9_withoutMap.txt",
	refgenome="mm9", 
	np=6, 
	mapdir='', # if mappability file directory is set to empty, mappability is not used. 
	atac=F     # Set TRUE for ATAC data. The default value is FALSE
	);

# This function call generates a hexamer bias frequency table from the given bamfile with mappability	
# This function call generates a hexamer bias frequency table from the given bamfile with mappability	
tabMappability = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_fed_mm9_withMap.txt", 
   refgenome="mm9", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM9);

