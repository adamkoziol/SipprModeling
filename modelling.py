__author__ = 'akoziol'

# Import the necessary modules
# OS is used for file/folder manipulations
import os
# Subprocess->call is used for making system calls
from subprocess import call
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob
# Shutil is useful for file moving/copying
import shutil
# prints variables in a form which can be used as input to the interpreter - similar to Data::Dumper?
import pprint
# Regex
import re

import time
import sys
import threading
import Queue
import shlex
import argparse



# Define the variables for the read length and fold coverage, respectively
#readLength = [30, 35, 40, 45, 50, 55, 60, 75, 80, 100, 150]
readLength = [54]
#foldCoverage = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100]
foldCoverage = [18]

# Initialize the required dictionaries
vcfData = {}

# Define the range of k-mer sizes for indexing of targets
#kmer = [5, 7, 9, 11, 13, 15]
kmer = [9]


os.chdir("/media/nas/akoziol/Pipeline_development/SipprModelling")
path = os.getcwd()
reference = "/media/nas/akoziol/Pipeline_development/SipprModelling/reference/Escherichia_coli_O157_H7_str_Sakai.fas"


def createSimulatedFiles(readLength, foldCoverage):
    """Iterates over the readLength and foldCoverage lists to create folders (if necessary)\
     and perform analyses"""
    os.chdir(path)
    for rLength in readLength:
        for fCov in foldCoverage:
            # Create a new folder(if necessary) at the appropriate location
            newPath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
            newFile = "%s/%s_%s_" % (newPath, rLength, fCov)
            artIlluminaCall = "art_illumina -i %s -l %s -f %s -m 225 -s 60 -o %s" % (reference, rLength, fCov, newFile)
            make_path(newPath)
            # Call art_illumina to simulate the reads into the appropriate folders
            # art_illumina -i /path-to-file/Escherichia_coli_O157_H7_str_Sakai.fas -l "readLength" -f "foldCoverage" \
            # -m 225 -s 60 -o /path-to-folder/Appropriate_name
            if not os.path.isfile("%s1.fq" % newFile):
                os.system(artIlluminaCall)


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def target_prep():
    """Prepares the target genes as necessary - this includes getting the files
    into a list, and processing using helper functions faidx_targets and index_targets"""
    os.chdir("%s/targets" % path)
    global targets
    targets = glob.glob("*.fa")
    faidx_targets()
    index_targets()


def index_targets():
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    for file in targets:
        # file.split splits the file (eg purA.fa) on the period. The first entry from the new list is stored in filename
        filename = file.split('.')[0]
        for size in kmer:
            # Create a new path to be created (if necessary) for the generation of the range of k-mers
            indexPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
            # Call the make_path function to make folders as necessary
            make_path(indexPath)
            os.chdir(indexPath)
            indexFile = "%s.smi" % filename
            if not os.path.isfile(indexFile):
                indexCommand = "smalt index -k %s -s 1 %s %s/targets/%s" % (size, filename, path, file)
                os.system(indexCommand)


def faidx_targets():
    """Creates .fai index files of the targets, which are necessary for the conversion
    of sorted BAM files to fastq files."""
    for file in targets:
        faidxFile = "%s.fai" % file
        faidxPath = "%s/targets/faidxFiles" % path
        make_path(faidxPath)
        if not os.path.isfile("%s/%s" % (faidxPath, faidxFile)):
            faidxCommand = "samtools faidx %s" % file
            os.system(faidxCommand)
            shutil.move(faidxFile, faidxPath)
            shutil.move(file, faidxPath)


def mapping(rLength, fCov, size, count, filename, megaName, filePath, total):
    """Performs the mapping of the simulated reads to the targets"""
    targetPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
    ###Include the commented out lines if necessary to have targets in the same directory as fastq files
    #os.chdir("%s/targets/%s/%s_%s" % (path, filename, filename, size))
    #for files in os.listdir("."):
    #    shutil.copy(files, newPath)
    if not os.path.isfile("%s.bam" % (megaName)):
        smaltMap = "smalt map -o %s.bam -n 24 -f bam -x -l pe %s/%s %s/%s_%s_1.fq %s/%s_%s_2.fq" \
                % (megaName, targetPath, filename, filePath, rLength, fCov, filePath, rLength, fCov)
        os.system(smaltMap)
        print "Mapping file %s of %s" % (count, total)


def sorting(count, megaName, sorted, sortedName, total):
    """Sorts the bam file in order for further manipulations to be possible"""
    # Sort the BAM file
    if not os.path.isfile(sortedName):
        bamSort = "samtools sort %s.bam %s" % (megaName, sorted)
        print "Sorting file %s of %s" % (count, total)
        os.system(bamSort)


def indexing(megaName, sortedName, count, total):
    """Indexes the sorted bam files in order to visualize the assemblies with tablet - note this is OPTIONAL"""
    indexedName = str(megaName) + "_sorted.bai"
    if not os.path.isfile(indexedName):
        bamIndex = "samtools index %s %s" % (sortedName, indexedName)
        print "Indexing file %s of %s" % (count, total)
        os.system(bamIndex)


def createVCF(target, sortedName, vcfFile, count, total):
    """Indexes the sorted bam files in order to visualize the assemblies with tablet - note this is OPTIONAL"""
    faidxTarget = "%s/targets/faidxFiles/%s" % (path, target)
    # Read this to understand why certain flags were used
    if not os.path.isfile(vcfFile):
        vcfPipe = "samtools mpileup -A -BQ0 -d 1000000 -uf %s %s | bcftools view -cg - > %s" % (faidxTarget, sortedName, vcfFile)
        print "Creating VCF file %s of %s" % (count, total)
        os.system(vcfPipe)


def populateVCFdata(megaName, vcfFile):
    """Opens the vcf file, extract the appropriate information, and places the data into the dictionary"""
    vcf = open(vcfFile, 'r')
        # Enter the first key in the dictionary - unlike Perl, Python doesn't support
    # autovivification  - intermediate data structures must be explicitly created
    vcfData[megaName] = {}
    for line in vcf:
        # All lines with a '#' are ignored
        if "#" in line:
            pass
        else:
            line = line.rstrip()
            tabLine = line.split('\t')
            pos = int(tabLine[1])
            ref = tabLine[3]
            alt = tabLine[4]
            qual = float(tabLine[5])
            info = tabLine[7]
            messyDP = info.split(';')[0]
            DP = messyDP.replace('DP=','')
            if alt == ".":
                base = ref
            else:
                base = alt
            # Populates the dictionary with the base position, base call, quality, and the depth of
            # coverage for that position the int() function had to be used here, as the sort function
            # used in traversing the dictionary was treating the values as strings, and sorting them
            # inappropriately eg (1, 10, 100, 11, 12, 13, ... 2) instead of (1, 2, ... 100)
            vcfData[megaName][pos] = {}
            vcfData[megaName][pos][base] = {}
            vcfData[megaName][pos][base][qual] = {}
            vcfData[megaName][pos][base][qual] = DP


def pipeline():
    """Calls the appropriate functions defined above and pipes them all together"""
    count = 0
    os.chdir(path)
    total = len(readLength) * len(foldCoverage) * len(targets) * len(kmer)
    # uses kmer, targets, readLength, foldCoverage
    for rLength in readLength:
        for fCov in foldCoverage:
            for target in targets:
                filename = target.split('.')[0]
                for size in kmer:
                    count += 1
                    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
                    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
                    newPath = "%s/%s" % (filePath, megaName)
                    sorted = str(megaName) + "_sorted"
                    sortedName = str(megaName) + "_sorted.bam"
                    os.chdir(newPath)
                    mapping(rLength, fCov, size, count, filename, megaName, filePath, total)
                    sorting(count, megaName, sorted, sortedName, total)
                    indexing(megaName, sortedName, count, total)
                    vcfFile = str(megaName) + "_sorted.vcf"
                    createVCF(target, sortedName, vcfFile, count, total)
                    outputFile = str(megaName) + "_output.csv"
                    populateVCFdata(megaName, vcfFile)


target_prep()
createSimulatedFiles(readLength, foldCoverage)
pipeline()
pprint.pprint(vcfData)

