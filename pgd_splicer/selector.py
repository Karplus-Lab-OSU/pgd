#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

from pgd.pgd_splicer.models import *

import urllib
import os
import re
import fileinput

class PDBSelectorTask():

    progressValue = 0

    def _work(self, args):
        if not os.path.exists(pdb_select_settings.PDB_TMP_DIR):
            os.mkdir(pdb_select_settings.PDB_TMP_DIR)

        file25 = pdb_select_settings.PDB_TMP_DIR+'/'+pdb_select_settings.PDB_SELECT_FILE_25
        file90 = pdb_select_settings.PDB_TMP_DIR+'/'+pdb_select_settings.PDB_SELECT_FILE_90

        #Download the files
        #getfile(pdb_select_settings.PDB_SELECT_URL+pdb_select_settings.PDB_SELECT_FILE_90, file90)
        #getfile(pdb_select_settings.PDB_SELECT_URL+pdb_select_settings.PDB_SELECT_FILE_25, file25)
        self.progressValue = 20

        #Parse select25 file
        select25 = self.parseFile(file25, pdb_select_settings.RESOLUTION)
        self.progressValue = 60

        #Merge select90 into select25 giving preference to select25
        merged = self.parseFile(file90, pdb_select_settings.RESOLUTION, select25)
        self.progressValue = 100

        #cleanup
        #os.remove(file25)
        #os.remove(file90)

        return merged

    """
    returns progress as a number between 0 and 100
    """
    def progress(self):
        return self.progressValue

    """
    Returns the status as a string
    """
    def progressMessage(self):
        return '%d of %d' % (self.count, self.stop)

    """
    Reset the task
    """
    def _reset(self):
        self.status = 0
        os.remove(file25)
        os.remove(file90)

    def parseFile(self, file, maxResolution, merge=None):
        if merge == None:
            map = {}
        else:
            map = merge

        #create regex pattern here so it is not done repeatedly while parsing file
        regexStr = '\s+(\d+)\s+([\w]+)\s+(\d+)\s+([-]?\d+\.\d+)\s+(\d+\.\d+)\s+([a-zA-Z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\w\-/(/) \[\]]+)'
        regexPattern = re.compile(regexStr)

        try:
            #iterate file
            for line in fileinput.input(file):

                if fileinput.lineno() < 4:
                    pass

                else:
                    protein = self.parsePDBSelectLine(line, regexPattern)
                    if not protein == None:
                        #Filter out proteins already in the list, this merges the lists giving preference to the filter
                        if not merge == None and protein[0] in merge:
                            continue

                        #check resolution range
                        if protein[2] > 0 and protein[2] <= maxResolution:
                            #print "Adding: %s" % protein
                            map[protein[0]] = protein


        finally:
            pass

        return map

    def parsePDBSelectLine(self, line, regexPattern):

        match = regexPattern.match(line)

        if match == None:
            return None

        '''protein.threshold = int(match.group(1))
        protein.identifier = match.group(2)
        #protein.chainId = (match.group(3))
        protein.aaCount = int(match.group(4))
        protein.resolution = float(match.group(5))
        protein.Rfactor = float(match.group(6))
        protein.method = match.group(7)
        protein.n_sid = int(match.group(8))
        protein.n_bck = int(match.group(9))
        protein.n_naa = int(match.group(10))
        protein.n_hlx = int(match.group(11))
        protein.n_bta = int(match.group(12))
        protein.compound = match.group(13)'''

        protein = [match.group(2)
                    ,int(match.group(1))
                    ,float(match.group(5))
                ]

        return protein

def getfile(url, output):
        #delete old version of file first
        if os.path.exists(output):
            os.remove(output)

        #get file
        file =  urllib.urlopen(url)
        raw = file.read()

        #write contents to file
        try:
            output = open(output, "w")
            output.write(raw)
        finally:
            output.close()

if __name__ == '__main__':
   task =  PDBSelectorTask('CommandLine PDBSelector')
   pdbs = task._work(None)
   print pdbs
