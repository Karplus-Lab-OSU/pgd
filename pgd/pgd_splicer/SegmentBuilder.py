
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

import math

from tasks.tasks import *
from pgd_core.models import *
from pgd_search.models import *

"""
Task that processes protein chains to create or update segment objects

This task requires that proteins, chains, and residues have already been created
"""
class SegmentBuilderTask(Task):

    def _work(self, args):
        length = 10
        lastIndex = length - 1

        #calculate the position of i in an array of size 'length'
        iIndex = int(math.ceil(length/2.0)-1)

        #calculate offset of last index from iIndex
        lastIndexOffset = lastIndex - iIndex

        proteins = Protein.objects.all()
        for protein in proteins:
            print 'protein: %s' % protein.code

            chains = protein.chains.all()
            for chain in chains:
                print '    chain: %s' % chain.code
 
                #setup initial list to have values where the first iteration of the list 
                #will be what is needed for processing the first residue
                segmentList = []

                #populate beginning of list with Nones
                for i in range(iIndex):
                    segmentList.append(None)

                #populate remaining residues with values looked up from the database
                id = 0
                for i in range(iIndex, length):
                    result = chain.residues.filter(chainIndex=id)
                    if len(result):
                        segmentList.append(result[0])
                    else:
                        segmentList.append(None)
                    id = id + 1

                print '        initial list: %s' % segmentList

                #determine index of the last known residue in the chain
                result = chain.residues.order_by('chainIndex').reverse()[0]
                chainLength = result.chainIndex
                print '        chainlength: %s' % chainLength

                #iterate through all possible residue indexes for this chain
                for ri in range(1, chainLength+1):

                    #roll list to next possible segment
                    #  this will roll lastIndexOffset past the last known index but this is expected
                    #  the queries will return None which is also expected
                    #  the last known index will have a maximum length of 1 with all remaining residues as None
                    result = chain.residues.filter(chainIndex=ri+lastIndexOffset)
                    if len(result):
                        segmentList.append(result[0])
                    else:
                        segmentList.append(None)
                    del segmentList[0]
                    #print '        list: %s' % segmentList

                    #if i is None then skip this segment
                    #its maxLength would be 0 and the segment would never be
                    #returned in any search
                    if not segmentList[iIndex]:
                        continue

                    #find existing segment or create new one
                    iProperty = 'r%d_index' % iIndex
                    segmentResult = Segment.objects.filter(protein__code = protein.code, i__chainIndex = segmentList[iIndex].chainIndex)
                    if len(segmentResult):
                        segment = segmentResult[0]
                    else:
                        segment = Segment()

                    #set residues in segment
                    for i in range(length):
                        pass
                        #segment.residue[i] = segmentList[i]
                    segment.protein = protein
                    segment.i = segmentList[iIndex]
 
                    #calculate max length of this particular segment
                    for i in range(iIndex):
                        # have we reached the max size for segments or found a None
                        if segmentList[iIndex+1] > lastIndex or segmentList[iIndex-1] < 0 or segmentList[iIndex+1] == None or segmentList[iIndex-1] == None:
                                segment.length = i-1
                                break

                    #save segment
                    #print '            segment: %s' % segment
                    #segment.save()

if __name__ == '__main__':
    builder = SegmentBuilderTask('Command Line Builder')
    builder._work(None)
