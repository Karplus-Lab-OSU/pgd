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


from ProcessPDBTask import *
from ftpupdate import *
from SegmentBuilder import *

class ProteinImportTask():
    """
    Full import process for a single protein.
    """
    description = 'Import task for a protein PDB.  Will download and process into both raw and searchable structures'

    def __init__(self, msg=None):
        TaskContainer.__init__(self, msg)

        self.add_task(FTPUpdateTask('Download required PDBs'))
        self.add_task(ProcessPDBTask('Process PDB into raw datastructure'))
        self.add_task(SegmentBuilderTask('Process raw structures into search table'))


class ParallelProteinImportTask():
    """
    Wrapper around ProteinImportTask to make it run in Parallel
    """
    description = 'Wrapper around ProteinImportTask to make it run in Parallel'

    def __init__(self):
        ParallelTask.__init__(self)
        self.subtask = ProteinImportTask()
        logger.debug( self.__dict__ )

    def _work(self, **kwargs):
        """
        overridden to batch up workunits
        """
        if kwargs and kwargs.has_key('data'):
            _data = kwargs['data']
            batch_size = int(kwargs['batch_size']) if kwargs.has_key('batch_size') else 50
            self.version = kwargs['version']

            logger.debug('ParallelProteinImportTask - repackaging work into batches of %s' % batch_size)
            kwargs['data'] = [_data[i:i+batch_size] for i in range(0, len(_data), batch_size)]
            logger.debug('ParallelProteinImportTask - %s workunits' % len(kwargs['data']))

        ParallelTask._work(self, **kwargs)

    def work_complete(self):
        set_setting_value('pgd_splicer.models','','DATA_VERSION',self.version)


if __name__ == '__main__':
    """
    When run from the command line this file should accept a list of protein
    codes and process them.

    Note: ProcessPDBTask currently excepts certain data to be included with PDB such as
          chains, threshold, resolution, rfactor, rfree.  Its not currently possible to
          pass those into this file.
    """
    import sys

    task = ProteinImportTask()

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
    task.work(**{'data':pdbs})
