#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

    # ==========================================================
    # Setup django environment 
    # ==========================================================
    if not os.environ.has_key('DJANGO_SETTINGS_MODULE'):
        os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
    # ==========================================================
    # Done setting up django environment
    # ==========================================================

from pydra_server.cluster.tasks import TaskContainer, ParallelTask
from ProcessPDBTask import *
from ftpupdate import *
from SegmentBuilder import *

class ProteinImportTask(TaskContainer):
    """
    Full import process for a single protein.
    """
    description = 'Import task for a protein PDB.  Will download and process into both raw and searchable structures'

    def __init__(self, msg=None):
        TaskContainer.__init__(self, msg)

        self.add_task(FTPUpdateTask('Download required PDBs'))
        self.add_task(ProcessPDBTask('Process PDB into raw datastructure'))
        self.add_task(SegmentBuilderTask('Process raw structures into search table'))


class ParallelProteinImportTask(ParallelTask):
    """
    Wrapper around ProteinImportTask to make it run in Parallel
    """
    description = 'Wrapper around ProteinImportTask to make it run in Parallel'

    def __init__(self, msg):
        ParallelTask.__init__(self, msg)
        self.subtask = ProteinImportTask()



if __name__ == '__main__':
    """
    When run from the command line this file should accept a list of protein
    codes and process them.

    Note: ProcessPDBTask currently excepts certain data to be included with PDB such as
          chains, threshold, resolution, rfactor, rfree.  Its not currently possible to
          pass those into this file.
    """
    from pydra_server.cluster.worker_proxy import WorkerProxy
    import sys

    task = ProteinImportTask()
    task.parent = WorkerProxy()

    pdbs = {}
    argv = sys.argv
    pdbs = []
    for i in range(1,len(argv),4):
        try:
            pdbs.append({'code':argv[i],
                      'threshold':float(argv[i+1]),
                      'resolution':float(argv[i+2]),
                      'rfactor':float(argv[i+3])
                      })
        except IndexError:
            print 'Usage: process_protein.py code threshold resolution rfactor ...'
            sys.exit(0)

    print pdbs
    task._work(**{'data':pdbs})