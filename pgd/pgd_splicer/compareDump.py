#!/usr/bin/env python

import fileinput
import sys
import math

def parseFile(file, codeIndex, residueIndex, LineFillChar):

    residues = []
    firstSegment = []

    #iterate php file building dict of values
    counting = 1
    for line in fileinput.input(file):

        if fileinput.lineno() < 2:
            pass

        #elif fileinput.lineno() > 10:
        #    break

        else:
            #first residue needs all rows processed so we can find the length of the segments
            if counting == 1:
                split = line.split('\t')

                #look for the first row in the second residue9
                if len(firstSegment) > 0 and firstSegment[-1][0] == LineFillChar and split[0] <> LineFillChar:
                    #we have all the residues in the first segment
                    #calculate length, iIndex and add the first residue    
                    segmentLength = len(firstSegment)
                    iIndex = int(math.ceil(segmentLength/2.0)-1)
                    iResidue = firstSegment[iIndex]                    
                    residues.append( '%s-%s' % (iResidue[codeIndex],iResidue[residueIndex]) ) 
                    print 'Segment Length: %s' % segmentLength

                    counting = 0;

                    #if the segment length is 1, then the row we read is 'i' otherwise it is throwaway data
                    if segmentLength == 1:
                        residues.append( '%s-%s' % (split[codeIndex],split[residueIndex]) )
                        count = 0
                    else:
                        count = 1
                    

                else:
                    # still on first residue keep adding data
                    #print 'append: %s' % split
                    #print len(split)
                    firstSegment.append(split)
            
            #all other residues we skip all rows except the row containing data for i
            else:
                if count == iIndex:
                    split = line.split('\t') 
                    residues.append( '%s-%s' % (split[codeIndex],split[residueIndex]) )   
                                
                count += 1

                #check for end of segment
                if count == segmentLength:
                    count = 0

    return segmentLength, residues
        
def compareDumps(php, python):
    php_length, php_residues = parseFile(php,1,3,' ')
    print 'PHP : Segment Length: %s' % php_length
    print '      Results Count : %s' % len(php_residues)
    print '-------------------------------'
    py_length, py_residues = parseFile(python,1,3,'')
    print 'PY  : Segment Length: %s' % py_length
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
    print '============================================='
    print 'Results missing from python'
    print '============================================='
    print not_in_py
    print    
    print '============================================='
    print 'Results that should not be in python results'
    print '============================================='
    print py_residues

if __name__ == '__main__':
    
    if len(sys.argv) < 3:
        print 'Usage: compareDump.py <php_dump> <python_dump>'
   
    compareDumps(sys.argv[1], sys.argv[2])
