Logic to the modeller:

Simulate Illumina reads using ART Illumina read simulator. Use Sakai as reference genome.
Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594 

Reads are simulated with variable lengths and fold coverage:
Flags used: 
-i (input file) - Sakai reference
-l (length of reads to be simulated)
-f (fold of read coverage to be simulated)
-m (mean length of reads to be simulated) - this was determined using pear_library_size_estimator.pl on the sequencing run from 2014-03-25
-s (standard deviation of reads to be simulated) - (see -m for how this was determined)
-o (output file)

art_illumina -i /path-to-file/Escherichia_coli_O157_H7_str_Sakai.fas -l "readLength" -f "foldCoverage" -m 225 -s 60 -o /path-to-folder/Appropriate_name

ReadLength
30
35
40
45
50
55
60
65
70
75
80
100
150
200
250

FoldCoverage
1
2
5
10
15
20
25
30
35
40
45
50
75
100

Targets to map will be full length vt1, vt2, eae, O157(O-Typer), adk(12), fumC(12), gyrB(8), icd(12), mdh(15), purA(2), recA(2) (MLST - appropriate alleles for Sakai).

Use GeneSippr approach to map raw reads to targets:
Index targets (smalt index, samtools faidx)
Reference Mapping (smalt map) - -x flag?
## ? ## -Use BAM output instead in smalt map - Replaces: Convert to BAM (samtools view)
Sort BAM file (samtools sort)
## ? ## samtools index -> samtools index rL50_fC20_fumC_kmer9_sorted.bam rL50_fC20_fumC_kmer9_sorted.bai
NOT necessary - replaced by VCF file mining. Use bedtools genomeCoverageBed to calculate the depth of coverage at each position of the target
	genomeCoverageBed -d -ibam rL35_fC50_fumC_kmer9_sorted.bam > output
Convert to VCF (custom pipe): samtools mpileup -uf fumC.fa rL50_fC10_fumC_kmer9_sorted.bam | bcftools view -cg - > rL50_fC10_fumC_kmer9_sorted.vcf
###? Covert to FASTQ: 
Convert to fasta (Biopython/python)
Find % identity (BLAST/python) - or some kind of variant calling?

Find maxima of accuracy/quality, time, price

Logic for fold coverage:
"Normal run" - OLC797 in run Ecoli_16_2013_08_28
12 072 757 reads PF, OLC797 accounts for 7.95% -> ~959 784 reads per direction -> ~ 2 M reads total

# 250 bp reads
foldCoverage = (numReads * readLength) / genomeSize
	     = 2 000 000 * 250 / 5 500 000
	     = 90.9

foldCoverage with rL of 50 = 18.2
foldCoverage with rL of 35 = 12.7

"Low run" - 2014-SEQ-173 in run 2014-04-11-Metagenome
10 205 375 reads PF, 2014-SEQ-173 accounts for 1.61% -> 164 306 reads per direction -> 328 613 reads total

# rL of 250
foldCoverage = 328 613 * 250 / 5 500 000
	     = 14.93

# rL of 50
foldCoverage = 328 613 * 50 / 5 500 000
	     = 2.99

"High run" - 2014-SEQ-070 in run 20140307-Metagenome
25 963 550 reads PF, 2014-SEQ-070 accounts for 11.13% -> 2 889 743 per direction -> 5 779 486 reads total

# rL of 250
foldCoverage = 5 779 486 * 250 / 5 500 000
	     = 262.7

#rL of 50
foldCoverage = 5 779 486 * 50 / 5 500 000
	     = 52.5

"Theoretical average"
Recent runs average PF = (17941796 + 10205375 + 17449780 + 17706900 + 16107341 + 16831848 + 25963550) / 7
		       = 122 206 590 / 7
		       = 17 458 084
standard deviation     = 4 609 002
standard error	       = 1 742 039	

"Theoretical average - no outliers"
Recent runs average PF = (17941796 + 17449780 + 20176796 + 17706900 + 16107341 + 16831848) / 6
		       = 86037665 / 5
		       = 17 702 410
standard deviation     = 1 381 728
standard error	       = 564 088	


