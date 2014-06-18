__author__ = 'blais'

import multiprocessing
import os
import sys
import subprocess
import errno
import glob
from multiprocessing import Pool

# Define the variables for the read length and fold coverage, respectively
readLength = [21, 22, 23, 24]
foldCoverage = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
kmer = [5, 7, 9]

os.chdir("/media/nas/akoziol/Pipeline_development/SipprModelling")
path = os.getcwd()
reference = "/media/nas/akoziol/Pipeline_development/SipprModelling/reference/Escherichia_coli_O157_H7_str_Sakai.fas"

os.chdir("%s/targets" % path)
targets = glob.glob("*.fa")


def createSimulatedFiles((rLength, fCov)):
    """Iterates over the readLength and foldCoverage lists to create folders (if necessary)\
     and perform analyses"""
    os.chdir(path)
    #for fCov in foldCoverage:
    # Create a new folder(if necessary) at the appropriate location
    newPath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newFile = "%s/%s_%s_" % (newPath, rLength, fCov)
    artIlluminaCall = "art_illumina -i %s -l %s -f %s -m 225 -s 60 -o %s" % (reference, rLength, fCov, newFile)
    make_path(newPath)
    # Call art_illumina to simulate the reads into the appropriate folders
    # art_illumina -i /path-to-file/Escherichia_coli_O157_H7_str_Sakai.fas -l "readLength" -f "foldCoverage" \
    # -m 225 -s 60 -o /path-to-folder/Appropriate_name
    if not os.path.isfile("%s1.fq" % newFile):
        #subprocess.Popen(shlex.split(artIlluminaCall), stdout=open(os.devnull, 'wb'))
        sys.stdout.write('.')
        subprocess.call(artIlluminaCall, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

    else:
        print sys.stdout.write('.')


def createSimulatedFilesProcesses(readLength, foldCoverage):
    args = []
    if __name__ == '__main__':
        pool = Pool()
        for rLength in readLength:
            for fCov in foldCoverage:
                args.append((rLength, fCov))
        pool.map(createSimulatedFiles, args)



def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

daemonProcess(readLength, foldCoverage)