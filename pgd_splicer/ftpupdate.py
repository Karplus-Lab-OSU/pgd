#!/usr/bin/env python

from __future__ import division

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

# from pgd_splicer.models import *
from pgd_splicer.models import ftp_update_settings

import os

from ftplib import FTP, error_perm
import re
import time

class FTPFile(object):

    def __init__(self, _flags, _wtf, _owner, _group, _size, _date, _name):
        self.flags = _flags
        self.wtf = _wtf
        self.owner = _owner
        self.group = _group
        self.size = _size
        self.date = _date
        self.name = _name

    def __str__(self):
        return '%s %s %s %s %s %s %s' % (self.flags, self.wtf, self.owner, self.group, self.size, self.date, self.name)


def processFile(str):
    p = re.compile( '([ldwrx\-]{9,11})\s+(\d+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\w+)\s+(\d+)\s+([0-9:]{4,5})\s+([a-zA-Z0-9_.\-]+)'  )
    m = p.match(str)

    if m==None:
        print "None"
        return

    # get pdb identifier
    pdb = m.group(9)[3:7]

    #check if we need this file
    if pdb not in pdb_local_files:
        return

    # extract date
    dateStr = m.group(6)+' '+m.group(7)+' '+m.group(8)

    # check for timestamp in place of year then parse date accordingly
    if dateStr.find(':') != -1:
        #no year in string so add it to the end
        date = time.strptime(dateStr + time.strftime(' %Y'), '%b %d %H:%M %Y' )
    else:
        date = time.strptime(dateStr, '%b %d %Y' )

    #file = FTPFile(m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), date, m.group(9))
    #files.append(file)
    # XXX how does it get over here?
    pdb_remote_files[pdb] = date


class FTPUpdateTask(object):

    pdbTotal = 0
    pdbCount = 0

    def work(self, data, **kwargs):
        print 'FTPUpdateTask - Starting' #, data, kwargs

        #create local directory if needed
        if not os.path.exists(ftp_update_settings.PDB_LOCAL_DIR):
            os.mkdir(ftp_update_settings.PDB_LOCAL_DIR)

        print '  FTP:', ftp_update_settings.PDB_FTP_HOST
        print '  REMOTE_DIR:', ftp_update_settings.PDB_REMOTE_DIR

        self.ftp = FTP(ftp_update_settings.PDB_FTP_HOST)
        #ftp.set_debuglevel(2) #set ftp debug level so all messages are shown
        self.ftp.login()
        self.ftp.cwd(ftp_update_settings.PDB_REMOTE_DIR)

        # accept data as either a list of proteins, or a single protein
        if isinstance(data, list):
            requested_pdbs = [d['code'][:4].lower() for d in data]
        else:
            requested_pdbs = [data['code'][:4].lower()]

        requested_pdbs.sort()

        self.pdbTotal = len(requested_pdbs)

        print 'Downloading %d PDB files to %s' % (self.pdbTotal,
                                                  ftp_update_settings.PDB_LOCAL_DIR)

        #get all the local timestamps
        for pdb in requested_pdbs:
            #===============================================
            #1 get list of local files matching the pdbs and the file modification date

            filename = 'pdb%s.ent.gz' % pdb
            path = os.path.join(ftp_update_settings.PDB_LOCAL_DIR, filename)
            if os.path.exists(path):
                date = time.gmtime(os.path.getmtime(path))
            else:
                date = None

            print "Local file:"
            if date:
                print '  - %s : %s' % (pdb, time.asctime(date))
            else:
                print '  - %s : None' % (pdb)


            # ===============================================
            #2 iterate through list and purge pdbs that are not new or newer

            print '    Checking Remote File:', filename

            try:
                resp = self.ftp.sendcmd('MDTM %s' % filename)
            except error_perm:
                # 550 error (permission error) results when a file does not
                # exist, remove pdb from list
                print 'File Not Found:', pdb, error_perm
                continue

            remote_date = time.strptime(resp[4:], '%Y%m%d%H%M%S')

            # keep only pdbs in the list that are new or newer, remove all others.
            if date and time.mktime(remote_date) <= time.mktime(date):
                # Already up to date; continue.
                continue

            # ===============================================
            #3 download files that are in the list
            self.download_pdb(pdb)

            self.pdbCount += 1

        self.ftp.quit()
        print "All finished; grabbed %d of %d (%d%%)" % (self.pdbCount,
                                                         self.pdbTotal,
                                                         self.progress())
        return kwargs

    def download_pdb(self, pdb):
        filename = 'pdb%s.ent.gz' % pdb
        local_filename = '%s/%s' % (ftp_update_settings.PDB_LOCAL_DIR,
                                    filename)

        # Grab size and pretty-print our progress.
        size = self.ftp.size(filename)
        if size:
            print "  %s - %s (%.2f KiB)" % (self.progressMessage(), pdb,
                                                size / 1024),
        else:
            print "  %s - %s" % (self.progressMessage(), pdb),

        #remove file if it exists already
        if os.path.exists(local_filename):
            os.remove(local_filename)

        self._incoming_file = open(local_filename,"w")
        self.ftp.retrbinary('RETR %s' % filename, self.processChunk)
        self._incoming_file.close()
        sys.stdout.write("\n")

    def processChunk(self, data):
        """
        Callback for ftp transfers.  This function will be called as chunks of
        data are received from the ftp server
        """
        sys.stdout.write(".")
        sys.stdout.flush()
        self._incoming_file.write(data)

    def progress(self):
        """
        returns progress as a number between 0 and 100
        """
        if self.pdbTotal == 0:
            return 100
        else:
            return int(self.pdbCount / self.pdbTotal * 100)

    def progressMessage(self):
        """
        Returns the status as a string
        """

        return "[%d/%d %d%%]" % (self.pdbCount, self.pdbTotal,
                                 self.progress())

    def _reset(self):
        """
        Reset the task
        """
        self.pdbCount = 0
        self.pdbTotal = 0


# code run if file executed from command line
if __name__ == '__main__':
    import sys

    pdbs = []
    argv = sys.argv
    if len(argv) == 1:
        print 'Usage:'
        print '   ftpupdate code [code ...]'
        print ''
        print '   <cmd> | ftpupdate --pipein'
        print '   piped protein codes must be separated by newlines'
        sys.exit(0)

    elif len(argv) == 2 and argv[1] == '--pipein':
        pdbs = [{'code':line} for line in sys.stdin]

    else:
        pdbs = [{'code':code} for code in sys.argv[1:]]

    task = FTPUpdateTask()
    task.work(data=pdbs)
