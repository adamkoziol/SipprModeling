#__author__ = 'blais'
__author__ = 'akoziol'

# Import the necessary modules
# OS is used for file/folder manipulations
import os
# Subprocess->call is used for making system calls
import subprocess
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob
# Shutil is useful for file moving/copying
import shutil
# prints variables in a form which can be used as input to the interpreter - similar to Data::Dumper?
#import pprint
# Regex
#import re

import time
import sys
from threading import Thread
from Queue import Queue
import shlex
#import argparse

# Define the variables for the read length and fold coverage, respectively
readLength = [30, 35, 40, 45, 50, 55, 60, 75, 80, 100, 150]
#readLength = [33]
foldCoverage = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100]
#foldCoverage = [1, 2, 3]

# Initialize the required dictionaries
vcfData = {}

# Define the range of k-mer sizes for indexing of targets
kmer = [5, 7, 9, 11, 13, 15]
#kmer = [5, 7, 9]


os.chdir("/media/nas/akoziol/Pipeline_development/SipprModelling/SE")
path = os.getcwd()
reference = "/media/nas/akoziol/Pipeline_development/SipprModelling/SE/reference/Escherichia_coli_O157_H7_str_Sakai.fas"


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

SimulatedFileQueue = Queue()


def createSimulatedFiles(SimulatedFileQueue):
    """Iterates over the readLength and foldCoverage lists to create folders (if necessary)\
     and perform analyses"""
    while True:
        os.chdir(path)
        rLength, fCov = SimulatedFileQueue.get()
        # Create a new folder(if necessary) at the appropriate location
        newPath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
        newFile = "%s/%s_%s_" % (newPath, rLength, fCov)
        # artIlluminaCall = "art_illumina -i %s -l %s -f %s -m 225 -s 60 -o %s" % (reference, rLength, fCov, newFile)
        artIlluminaCall = "art_illumina -i %s -l %s -f %s -o %s" % (reference, rLength, fCov, newFile)
        make_path(newPath)
        # Call art_illumina to simulate the reads into the appropriate folders
        # art_illumina -i /path-to-file/Escherichia_coli_O157_H7_str_Sakai.fas -l "readLength" -f "foldCoverage" \
        # -m 225 -s 60 -o /path-to-folder/Appropriate_name
        if not os.path.isfile("%s.fq" % newFile):
            #subprocess.Popen(shlex.split(artIlluminaCall), stdout=open(os.devnull, 'wb'))
            #subprocess.call(artIlluminaCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('wub')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        SimulatedFileQueue.task_done()


def createSimulatedFilesThreads():
    """ ghkwjhg """
    print "Creating simulated files"
    for i in range(len(readLength) * len(foldCoverage)):
        threads = Thread(target=createSimulatedFiles, args=(SimulatedFileQueue,))
        threads.setDaemon(True)
        threads.start()
    for rLength in readLength:
        for fCov in foldCoverage:
            SimulatedFileQueue.put((rLength, fCov))

    #wait on the queue until everything has been processed
    SimulatedFileQueue.join()

faidxQueue = Queue()


def faidxTargets(faidxQueue):
    """Creates .fai index files of the targets, which are necessary for the conversion
    of sorted BAM files to fastq files."""
    while True:
        file = faidxQueue.get()
        faidxFile = "%s.fai" % file
        faidxPath = "%s/targets/faidxFiles" % path
        make_path(faidxPath)
        if not os.path.isfile("%s/%s" % (faidxPath, faidxFile)):
            faidxCommand = "samtools faidx %s" % file
            #subprocess.Popen(shlex.split(faidxCommand), stdout=open(os.devnull, 'wb'))
            subprocess.call(faidxCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            shutil.move(faidxFile, faidxPath)
            shutil.move(file, faidxPath)
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        faidxQueue.task_done()


def faidxTargetsThreads(targets):
    print '\nProcessing targets with faidx'
    for i in range(len(targets)):
        threads = Thread(target=faidxTargets, args=(faidxQueue,))
        threads.setDaemon(True)
        threads.start()
    for target in targets:
        faidxQueue.put(target)

    #wait on the queue until everything has been processed
    faidxQueue.join()

indexQueue = Queue()


def indexTargets(indexQueue):
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    while True:
        target, size = indexQueue.get()
        filename = target.split('.')[0]
        # Create a new path to be created (if necessary) for the generation of the range of k-mers
        indexPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
        # Call the make_path function to make folders as necessary
        make_path(indexPath)
        indexFileSMI = "%s.smi" % filename
        indexFileSMA = "%s.sma" % filename
        if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
            indexCommand = "smalt index -k %s -s 1 %s %s/targets/%s" % (size, filename, path, target)
            #subprocess.Popen(shlex.split(indexCommand), stdout=open(os.devnull, 'wb'))
            subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            shutil.move(indexFileSMI, indexPath)
            shutil.move(indexFileSMA, indexPath)
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        indexQueue.task_done()


def indexTargetsThreads(targets):
    #print ''
    print '\nIndexing targets'
    for i in range(len(targets)):
        threads = Thread(target=indexTargets, args=(indexQueue,))
        threads.setDaemon(True)
        threads.start()
    for target in targets:
        for size in kmer:
            indexQueue.put((target, size))

    #wait on the queue until everything has been processed
    indexQueue.join()


mappingQueue = Queue()

def mapping(mappingQueue):
    """Performs the mapping of the simulated reads to the targets"""
    while True:
        targets, rLength, fCov, target, size = mappingQueue.get()
        filename = target.split('.')[0]
        megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
        filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
        newPath = "%s/%s" % (filePath, megaName)
        make_path(newPath)
        targetPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
        if not os.path.isfile("%s/%s.bam" % (newPath, megaName)):
            # smaltMap = "smalt map -o %s/%s.bam -n 24 -f bam -x -l pe %s/%s %s/%s_%s_1.fq %s/%s_%s_2.fq" \
            #            % (newPath, megaName, targetPath, filename, filePath, rLength, fCov, filePath, rLength, fCov)
            smaltMap = "smalt map -o %s/%s.bam -f bam -x %s/%s %s/%s_%s_.fq" \
                       % (newPath, megaName, targetPath, filename, filePath, rLength, fCov)
            #subprocess.Popen(shlex.split(smaltMap), stdout=open(os.devnull, 'wb'))
            subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        mappingQueue.task_done()


def mappingThreads(targets):
    """Mapping threads!"""
    os.chdir(path)
    print '\nPerforming reference mapping'
    length = len(readLength) * len(foldCoverage * len(targets) * len(kmer))
    for i in range(length):
        threads = Thread(target=mapping, args=(mappingQueue,))
        threads.setDaemon(True)
        threads.start()
    # uses kmer, targets, readLength, foldCoverage
    for rLength in readLength:
        for fCov in foldCoverage:
            for target in targets:
                for size in kmer:
                    mappingQueue.put((targets, rLength, fCov, target, size))
            mappingQueue.join()
    #wait on the queue until everything has been processed
    mappingQueue.join()

sortingQueue = Queue()


def sorting(sortingQueue):
    """Sorts the bam file in order for further manipulations to be possible"""
    while True:
        targets, rLength, fCov, target, size = sortingQueue.get()
        filename = target.split('.')[0]
        megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
        sorted = megaName + "_sorted"
        sortedMegaName = megaName + "_sorted.bam"
        filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
        newPath = "%s/%s" % (filePath, megaName)
        make_path(newPath)
        #Sort the BAM file
        if not os.path.isfile("%s/%s" % (newPath, sortedMegaName)):
            bamSort = "samtools sort %s/%s.bam %s/%s" % (newPath, megaName, newPath, sorted)
            #subprocess.Popen(shlex.split(bamSort), stdout=open(os.devnull, 'wb'))
            subprocess.call(bamSort, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')

        # signals to queue job is done
        sortingQueue.task_done()

def sortingThreads(targets):
    print '\nSorting bam files'
    length = len(readLength) * len(foldCoverage * len(targets) * len(kmer))
    for i in range(length):
        threads = Thread(target=sorting, args=(sortingQueue,))
        #threads = Sorting(sortingQueue)
        threads.setDaemon(True)
        threads.start()
        # uses kmer, targets, readLength, foldCoverage
    for rLength in readLength:
        for fCov in foldCoverage:
            for target in targets:
                for size in kmer:
                    sortingQueue.put((targets, rLength, fCov, target, size))

            sortingQueue.join()
    #wait on the queue until everything has been processed
    sortingQueue.join()

bamIndexQueue = Queue()


def bamIndexing(bamIndexQueue):
    """Indexes the sorted bam files in order to visualize the assemblies with tablet - note this is OPTIONAL"""
    while True:
        targets, rLength, fCov, target, size = bamIndexQueue.get()
        filename = target.split('.')[0]
        megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
        sortedMegaName = megaName + "_sorted.bam"
        filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
        newPath = "%s/%s" % (filePath, megaName)
        indexedName = megaName + "_sorted.bai"
        if not os.path.isfile("%s/%s" % (newPath, indexedName)):
            bamIndex = "samtools index %s/%s %s/%s" % (newPath, sortedMegaName, newPath, indexedName)
            #subprocess.Popen(shlex.split(bamIndex), stdout=open(os.devnull, 'wb'))
            subprocess.call(bamIndex, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        bamIndexQueue.task_done()


def bamIndexingThreads(targets):
    print '\nIndexing bam files'
    length = len(readLength) * len(foldCoverage * len(targets) * len(kmer))
    for i in range(length):
        threads = Thread(target=bamIndexing, args=(bamIndexQueue,))
        threads.setDaemon(True)
        threads.start()
        # uses kmer, targets, readLength, foldCoverage
    for rLength in readLength:
        for fCov in foldCoverage:
            for target in targets:
                for size in kmer:
                    bamIndexQueue.put((targets, rLength, fCov, target, size))

    #wait on the queue until everything has been processed
    bamIndexQueue.join()

vcfQueue = Queue()


def createVCF(vcfQueue):
    """"""
    while True:
        rLength, fCov, target, size = vcfQueue.get()
        filename = target.split('.')[0]
        megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
        sortedMegaName = megaName + "_sorted.bam"
        filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
        vcfFile = megaName + "_sorted.vcf"
        newPath = "%s/%s" % (filePath, megaName)
        faidxTarget = "%s/targets/faidxFiles/%s" % (path, target)
        # Read this to understand why certain flags were used
        # http://samtools.sourceforge.net/mpileup.shtml
        if not os.path.isfile("%s/%s" % (newPath, vcfFile)):
            vcfPipe = "samtools mpileup -A -BQ0 -d 1000000 -uf %s %s/%s | bcftools view -cg - > %s/%s" \
                      % (faidxTarget, newPath, sortedMegaName, newPath, vcfFile)
            #print vcfPipe
            #os.system(vcfPipe)
            #subprocess.Popen(shlex.split(vcfPipe), stdout=open(os.devnull, 'wb'))
            subprocess.call(vcfPipe, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
            sys.stdout.write('.')
        else:
            sys.stdout.write('.')
        # signals to queue job is done
        vcfQueue.task_done()


def createVCFThreads(targets):
    print '\nCreating vcf files'
    length = len(readLength) * len(foldCoverage * len(targets) * len(kmer))
    for i in range(length):
        threads = Thread(target=createVCF, args=(vcfQueue,))
        threads.setDaemon(True)
        threads.start()
        # uses kmer, targets, readLength, foldCoverage
    for rLength in readLength:
        for fCov in foldCoverage:
            for target in targets:
                for size in kmer:
                    vcfQueue.put((rLength, fCov, target, size))
            vcfQueue.join()
    #wait on the queue until everything has been processed
    vcfQueue.join()

def pipeline():
    """Calls all the functions in a way that they can be multithreaded"""
    createSimulatedFilesThreads()
    os.chdir("%s/targets" % path)
    targets = glob.glob("*.fa")
    faidxTargetsThreads(targets)
    indexTargetsThreads(targets)
    # Start the mapping operations
    mappingThreads(targets)
    sortingThreads(targets)
    bamIndexingThreads(targets)
    createVCFThreads(targets)
    #

start = time.time()
pipeline()
print "\nElapsed Time: %s seconds" % (time.time() - start)