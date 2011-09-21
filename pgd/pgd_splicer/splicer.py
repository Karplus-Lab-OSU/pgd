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

from django import forms

from dunbrack_selector import DunbrackPDBSelectorTask
from process_protein import ParallelProteinImportTask


class SplicerInput(forms.Form):
    """
    Form object used to render an interface that captures input
    for TestTask.
    """
    thresholds = forms.CharField(initial='25,90', help_text='Enter a comma delimeted list of thresholds to include')
    resolution = forms.FloatField(label='Max Resolution', initial='3.0', help_text='Enter the max resolution to accept')
    rfactor    = forms.FloatField(label='Max R-Factor', initial='1.0', help_text='Enter the max r-factor to accept')
    batch_size = forms.IntegerField(label='Batch Size', initial='50', help_text='Enter the number of proteins to include per workunit.  A protein takes between 1 and 15 seconds to process depending on chain count.  Proteins that are up to date will take only about 1 second or less to process.  Batches should ideally be a few minutes long.  A good benchmark is 25-50 for a full import.  1000 for a weekly update with only a few hundred proteins to update.')


class SplicerTask(TaskContainer):
    """
    Splicer - The import tool for the Protein Geometry Database (PGD)

    Selects, Downloads, and Imports protein geometry from PDB files.
    """
    description = 'Splicer - The import tool for the Protein Geometry Database (PGD)'
    form = SplicerInput

    def __init__(self, msg=None):
        TaskContainer.__init__(self, msg)

        self.add_task(DunbrackPDBSelectorTask('Select PDBs for download'),1)
        self.add_task(ParallelProteinImportTask())