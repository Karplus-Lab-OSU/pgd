import urllib
import os
import re
import fileinput


RESOLUTION=1.75                                                 # Resolution cutoff value
PDB_DIR='pdb                                                     # directory of pdb files'
PDB_SELECT_URL='http://bioinfo.tg.fh-giessen.de/pdbselect/'     #URL for pdbselect
PDB_SELECT_URL_FILE_90='recent.pdb_select90'                    # file to retreive
PDB_SELECT_URL_FILE_25='recent.pdb_select25'                    # file to retreive
PDB_TMP_DIR='tmp'
PDB_SELECT_FILE_90='recent.pdb_select90'                   # Filename for pdbselect
PDB_SELECT_FILE_25='recent.pdb_select25'


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


def createDownloadList():
    if not os.path.exists(PDB_TMP_DIR):
        os.mkdir(PDB_TMP_DIR)

    file25 = PDB_TMP_DIR+'/'+PDB_SELECT_FILE_25
    file90 = PDB_TMP_DIR+'/'+PDB_SELECT_FILE_90

    #Download the files
    #getfile(PDB_SELECT_URL+PDB_SELECT_URL_FILE_90, file90)
    #getfile(PDB_SELECT_URL+PDB_SELECT_URL_FILE_90, file25)

    #Parse select25 file
    
    select25 = parseFile(file25, RESOLUTION)
    print len(select25)

    #Merge select90 into select25 giving preference to select25
    merged = parseFile(file90, RESOLUTION, select25)

    #cleanup
    #os.remove(file25)
    #os.remove(file90)

    return merged


def parseFile(file, maxResolution, merge=None):
    if merge == None:
        map = {}
    else:
        map = merge

    #create regex pattern here so it is not done repeatedly while parsing file
    regexStr = '\s+(\d+)\s+([\w]+)\s+(\d+)\s+([-]?\d+\.\d+)\s+(\d+\.\d+)\s+([a-zA-Z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\w\-/(/) \[\]]+)'
    regexPattern = re.compile(regexStr)

    try:
        #open file
        #file = open(file, "r")

        #iterate file
        for line in fileinput.input(file):

            if fileinput.lineno() < 4:
                pass

            else:
                protein = parsePDBSelectLine(line, regexPattern)
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
        #file.close()

    return map

def parsePDBSelectLine(line, regexPattern):

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

if __name__ == '__main__':
   pdbs =  createDownloadList()
   print len(pdbs)
