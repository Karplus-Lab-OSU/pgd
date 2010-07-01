#!/usr/bin/env python

import fileinput
import sys
import math

def parseFile(file, codeIndex, residueIndex, php=False):

    residues = []
    firstSegment = []

    #iterate php file building dict of values
    first = True
    meta_end = 1000
    for line in fileinput.input(file):

        if fileinput.lineno() < meta_end:
            if line == '***END_META_DATA***':
                meta_end = fileinput.lineno()+1
            pass

        #elif fileinput.lineno() > 10:
        #    break

        else:
            #first residue needs all rows processed so we can find the length of the segments
            if first:
                split = line.split('\t')

                #look for the first row in the second residue9
                if len(firstSegment) > 0 and split[0] <> firstSegment[0][0]:
                    
                    #we have all the residues in the first segment
                    #calculate length, iIndex and add the first residue
                    segmentLength = len(firstSegment)
                    iIndex = int(math.ceil(segmentLength/2.0)-1)
                    iResidue = firstSegment[iIndex]
                    
                    residues.append( '%s-%s' % (iResidue[codeIndex],iResidue[residueIndex]) ) 

                    first = False;

                    #if the segment length is 1, then the row we read is 'i' otherwise it is throwaway data
                    if segmentLength == 1:
                        residues.append( '%s-%s' % (split[codeIndex],split[residueIndex]))
                        count = 0
                    else:
                        count = 1


                else:
                    # still on first residue keep adding data
                    firstSegment.append(split)


            #all other residues we skip all rows except the row containing data for i
            else:
                if count == iIndex:
                    split = line.split('\t') 
                    residues.append( '%s-%s' % (split[codeIndex],split[residueIndex]))

                count += 1

                #check for end of segment
                if count == segmentLength:
                    count = 0

    return segmentLength, residues

def compareDumps(php, python):
    php_length, php_residues = parseFile(php,1,3)
    print 'Old : Segment Length: %s' % php_length
    print '      Results Count : %s' % len(php_residues)
    print '-------------------------------'
    py_length, py_residues = parseFile(python,1,3)
    print 'New  : Segment Length: %s' % py_length
    print '      Results Count : %s' % len(py_residues)
    print '-------------------------------'
    not_in_py = []

    for code in php_residues:
        if code in py_residues:
            py_residues.remove(code)

        else:
            not_in_py.append(code)

    print 'Residues missing from python results: %s' % len(not_in_py)
    print 'Residues that should not be in python results: %s' % len(py_residues)

    ignore = [
        '3DUW','1KNM','1WM3','2RFR',
        '1IXH','1O5X','2FE5','2IIM','2OLN','2VVG','3E8M'
    ]


    if len(not_in_py):
        print '============================================='
        print 'Results missing from new Dump'
        print '============================================='
        for i in not_in_py:
            if not i[:4] in ignore:
                print i

    return 
    if len(py_residues):
        print
        print '============================================='
        print 'Results that should not be in new dump results'
        print '============================================='
        print py_residues

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print 'Usage: compareDump.py <old_dump> <new_dump>'

    compareDumps(sys.argv[1], sys.argv[2])
