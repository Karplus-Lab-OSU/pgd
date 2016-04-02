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
        os.environ['DJANGO_SETTINGS_MODULE'] = 'pgd.settings'
    # ==========================================================
    # Done setting up django environment
    # ==========================================================


import gzip
import os
import re
import sys
import urllib

from django.conf import settings

def printc(txt):
    """ direct print statement that can be overritten to silence output """
    sys.stdout.write(txt)
    sys.stdout.write('\n')


class DunbrackPDBSelectorTask():
    """
    DunbrackPDBSelectorTask downloads and processes PDB selection files from
    Dunbrack's PISCE service.  It pulls files from the list of pre-culled files
    rather than performing a custom query.  After downloading the files, they
    are parsed and merged to form a single dictionary of PDBs and intitial
    data.
    """

    progressValue = 0

    def work(self, **kwargs):
        """
        Main work function.
        """

        printc("DunbrackPDBSelectorTask: Starting")
        try:
            thresholds = kwargs['thresholds']
            if isinstance(thresholds, (str,unicode)):
                thresholds = [int(s) for s in thresholds.split(',')]
            elif isinstance(thresholds, int):
                thresholds = [thresholds]

        except KeyError:
            thresholds = [25,90]

        try:
            resolution = kwargs['resolution']
        except KeyError:
            resolution = '3.0'

        try:
            r_factor = kwargs['r_factor']
        except KeyError:
            r_factor = '1.0'

        printc('  thresholds:     %s' % thresholds)
        printc('  max resolution: %s' % resolution)
        printc('  r_factor:       %s' % r_factor)

        printc('Saving pdbs in %s' % settings.PDB_TMP_DIR)
        if not os.path.exists(settings.PDB_TMP_DIR):
            os.mkdir(settings.PDB_TMP_DIR)
            printc('     making dir %s' % settings.PDB_TMP_DIR)

        # Download the files
        dunbrack_url = 'http://dunbrack.fccc.edu/Guoli/culledpdb_hh/'
        files = get_files(dunbrack_url, thresholds, resolution, r_factor)

        printc('files: %s %s %s %s' % (dunbrack_url, thresholds, resolution, r_factor))

        self.progressValue = 20

        # Parse the files that were downloaded - the downloaded files
        # might not match expected files if nothing matched the parameters
        # given to get_files
        proteins = {}
        for filename, threshold in files:
            path = os.path.join(settings.PDB_TMP_DIR, filename)
            proteins = self.parse_file(path, float(resolution), threshold, proteins)
            os.remove(path)

        # store the date from the first file to use as the version.  The version will be
        # updated now even though the import has just begun.  Its marked to indicate that
        # it is still in progress
        date = files[0][0][26:32]
        version = '20%s-%s-%s' % (date[:2], date[2:4], date[4:])
        #set_setting_value('pgd_splicer.models','','DATA_VERSION','%s (import in progress)' % version)

        self.progressValue = 100

        printc('VERSION: %s' % version)
        return {'version':version, 'data':[v for k,v in proteins.items()]}


    def progress(self):
        """
        returns progress as a number between 0 and 100
        """
        return self.progressValue


    def _reset(self):
        """
        Reset the task
        """
        self.status = 0


    def parse_file(self, path, max_resolution, threshold, proteins={}):
        """
        Parse a file from Dunbrack's PICSE server.  This will
        iterate all lines in the file and extract the list of
        proteins that match the resolution criteria.
        """

        """
            create regex pattern here so it is not done repeatedly while parsing file

            groups:
                0 - Protein ID
                1 - Chain ID
                2 - Length of protein chain
                3 - Exptl.
                4 - Resolution
                5 - R-factor
                6 - FreeRValue
        """
        regex_str = '(\w{4})(\w)\s+(\d+)\s+(\w+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)'
        regex_pattern = re.compile(regex_str)

        printc('Processing: %s' % path)

        raw = None
        try:
            _file = gzip.open(path, 'r')

            #first line is labels, discard it
            _file.readline()

            for line in _file:
                match = regex_pattern.match(line)
                if match:
                    groups = match.groups()

                    if groups[0] in proteins:
                        # if protein already exists just update the additional
                        # chain information.  The properties should not change
                        # between records in the selection file.
                        protein = proteins[groups[0]]
                        if not groups[1] in protein['chains']:
                            protein['chains'].append(groups[1])
                            #print 'Selecting Protein: %s   Chain: %s   Threshold: %s' % (groups[0],groups[1], threshold)

                    else:
                        # protein is not in proteins dict yet create initial
                        # structure from parsed properties.
                        resolution = float(groups[4])
                        if resolution > 0 and resolution <= max_resolution:
                            proteins[groups[0]] = {
                                'code':groups[0],
                                'chains':[groups[1]],
                                'resolution':groups[4],
                                'rfactor':groups[5],
                                'rfree':groups[6],
                                'threshold':threshold
                            }

                        #print 'Selecting Protein: %s   Chain: %s   Threshold: %s' % (groups[0],groups[1], threshold)

        finally:
            if _file:
                _file.close()

        return proteins


def get_files(url, thresholds, resolution, r_factor):
    """
    Download files from the given url that match the input properties.
    """
    selection_page = urllib.urlopen(url).read()

    threshold = '|'.join(['%s'% t for t in thresholds])
    pattern = '\s(cullpdb_pc(%s)_res%s_R%s_.*\d\.gz)' % (threshold, resolution, r_factor)

    files = re.findall(pattern, selection_page)

    output = None
    for filename,threshold in files:
        printc('Downloading: %s' % filename)
        #get file
        file =  urllib.urlopen(url +'/'+ filename )
        raw = file.read()

        #write contents to file
        try:
            output = open(os.path.join(settings.PDB_TMP_DIR, filename), "w")
            output.write(raw)

        finally:
            if output:
                output.close()

    return files


if __name__ == '__main__':
    import sys

    task =  DunbrackPDBSelectorTask()
    pipe = False

    argv = sys.argv

    if '--help' in argv:
        print 'Usage:'
        print '   dunbrack_selector.py [--file name] [--pipeout]'
        print ''
        print '   Options:'
        print '      --help - this screen'
        print '      --pipeout - output will skip header information so that it may be piped easier into other tools'
        print '      --codes - output only protein codes (default includes chains, threshold, resolution, rfactor, rfree)'
        sys.exit(0)
    elif '--pipeout' in argv:
        def null_print(txt):
            pass
        printc = null_print

    pdbs = task.work(clean=pipe)

    if '--codes' in argv:
        format = '%(code)s\n'
    else:
        format = '%(code)s %(chains)s %(threshold)s %(resolution)s %(rfactor)s %(rfree)s\n'

    out = sys.stdout
    for p in pdbs['data']:
        p['chains'] = ''.join(p['chains'])
        out.write(format % p)
