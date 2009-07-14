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

from pydra_server.cluster.tasks import TaskContainer
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


class SplicerTask(TaskContainer):
    """
    Splicer - The import tool for the Protein Geometry Database (PGD)

    Selects, Downloads, and Imports protein geometry from PDB files.
    """
    description = 'Splicer - The import tool for the Protein Geometry Database (PGD)'
    form = SplicerInput

    def __init__(self, msg=None):
        TaskContainer.__init__(self, msg)

        self.add_task(DunbrackPDBSelectorTask('Select PDBs for download'))
        self.add_task(ParallelProteinImportTask('ProcessSelectedProteins'))