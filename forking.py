from multiprocessing import Pool
import os
import sys
import subprocess
import errno
import glob

# def f(x):
#     return x*x
#
# if __name__ == '__main__':
#     pool = Pool(processes=4)              # start 4 worker processes
#     result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
#     print result.get(timeout=1)           # prints "100" unless your computer is *very* slow
#     print pool.map(f, range(10))          # prints "[0, 1, 4,..., 81]"

# Define the variables for the read length and fold coverage, respectively
readLength = [33]
foldCoverage = [10, 12]
kmer = [5, 7, 9]

os.chdir("/media/nas/akoziol/Pipeline_development/SipprModelling")
path = os.getcwd()
reference = "/media/nas/akoziol/Pipeline_development/SipprModelling/reference/Escherichia_coli_O157_H7_str_Sakai.fas"

os.chdir("%s/targets" % path)
targets = glob.glob("*.fa")


def createSimulatedFiles(rLength):
    """Iterates over the readLength and foldCoverage lists to create folders (if necessary)\
     and perform analyses"""

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
            #subprocess.Popen(shlex.split(artIlluminaCall), stdout=open(os.devnull, 'wb'))
            sys.stdout.write('.')
            #subprocess.call(artIlluminaCall, shell=True, stdout=open(os.devnull, 'wb'))
        else:
            print sys.stdout.write('.')


os.chdir(path)
if __name__ == '__main__':
    with Pool(processes=4) as pool:
        #pool.apply_async(createSimulatedFiles, readLength)
        pool.map(createSimulatedFiles, readLength)
        #pool.start()
        #pool.join()

createSimulatedFilesThreads()

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise