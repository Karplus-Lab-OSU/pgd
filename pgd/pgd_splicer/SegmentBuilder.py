from tasks.tasks import *
from pgd_core.models import *

"""
Task that processes protein chains to create or update segment objects

This task requires that proteins, chains, and residues have already been created
"""
class SegmentBuilderTask(Task):

    def getChainIds(self, code):
        def my_custom_sql(self):
        from django.db import connection
        cursor = connection.cursor()
        cursor.execute("SELECT distinct chainID FROM pgd_core_residue WHERE code_id = %s", [code])

        ids = []
        for row in cursor:
            ids.append(row)

        return ids

    def _work(self):
        length = 10
        lastIndex = length - 1

        #calculate the position of i in an array of size 'length'
        iIndex = int(math.ceil(length/2.0)-1)

        proteins = Protein.objects.all()
        for protein in protein:

            chains = protein.chains.all()
            for chain in chains:            
 
                #setup initial list to have values where the first iteration of the list 
                #will be what is needed for processing the first residue
                segmentList = []

                #populate beginning of list with Nones
                for i in range(iIndex+1):
                    segmentList[i] = None
                
                #populate remaining residues with values looked up from the database
                id = 0
                for i in range(iIndex+1:length):
                    result = chain.residues.filter(newId=id)
                    if len(result):
                        segmentList[i] = result[0]
                    else:
                        segmentList[i] = None
                    id++

                #determine max length of chain.  
                chainLength = ???

                #iterate through all possible residue id's for this chain
                for ri in range(id-1, chainLength):
                    
                    #roll list to next possible segment
                    result = chain.residues.filter(newId=ri)
                    if len(result):
                        segmentList.append(result[0])
                    else:
                        segmentList.append(None)
                       
                    segmentList[0]
                    id++

                    #if i is None then skip this segment
                    if not segmentList[iIndex]:
                        continue
                 
                    #find existing segment or create new one
                    segmentResult = Segment.objects.filter(iProperty = segmentList[iIndex].id, protein_id = protein.id)
                    if len(segmentResult):
                        segment = segmentResult[0]
                    else:
                        segment = Segment()

                    #set values
                    for i in range(length):
                        segment.residue[i] = segmentList[i]
                    segment.protein = protein
 
                    #calculate max length
                    for i in range():
                        
                        # have we reached the max size for segments or found any Nones
                        if segmentList[iIndex+1] > lastIndex or segmentList[iIndex-1] < 0 
                            or segmentList[iIndex+1 == None or segmentList[iIndex-1] == None:
                                segment.length = i-1
                                break

                    #save segment
                    segment.save()
