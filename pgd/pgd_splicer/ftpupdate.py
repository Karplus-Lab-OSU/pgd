#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    import os

    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

from tasks.tasks import *
from pgd_splicer.models import *

import os

from sets import Set
from ftplib import FTP, error_perm
import re
import time

class FTPFile():
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
    pdb_remote_files[pdb] = date


class FTPUpdateTask(Task):

    pdbTotal = 0
    pdbCount = 0
    processing_dates = None

    def _work(self, args):
        pdb_local_files = {}
        pdb_remote_files = {}

        #===============================================
        #1 get list of local files matching the pdbs and the file modification date
        self.processing_dates = True

        # build list of pdbs from keys
        requested_pdbs = []
        for key in args.keys():
            #strip chain indicator
            requested_pdbs.append(key[:4].lower())

        #create local directory if needed
        if not os.path.exists(ftp_update_settings.PDB_LOCAL_DIR):
            os.mkdir(ftp_update_settings.PDB_LOCAL_DIR)

        #get all the local timestamps
        for pdb in requested_pdbs:
            path = '%s/pdb%s.ent.Z' % (ftp_update_settings.PDB_LOCAL_DIR, pdb)
            if os.path.exists(path):
                date = time.gmtime(os.path.getmtime(path))
            else:
                date = None

            pdb_local_files[pdb.lower()] = date


        print "LOCAL FILES:"
        for entry in pdb_local_files:
            if pdb_local_files[entry] == None:
                print '  - %s : None' % (entry)
            else:
                print '  - %s : %s' % (entry, time.asctime(pdb_local_files[entry]))

        print "------------------------------------------"

        # ===============================================
        #2 iterate through list and purge pdbs that are not new or newer
        ftp = FTP(ftp_update_settings.PDB_FTP_HOST)
        #ftp.set_debuglevel(2) #set ftp debug level so all messages are shown
        ftp.login()
        ftp.cwd(ftp_update_settings.PDB_REMOTE_DIR)

        for pdb in pdb_local_files:
            filename = 'pdb%s.ent.Z' % pdb

            localdate = pdb_local_files[pdb]
            try:
                resp = ftp.sendcmd('MDTM %s' % filename)
                #date = time.strptime('%s%s' % (resp[4:],time.tzname[0]), '%Y%m%d%H%M%S%Z' )
                date = time.strptime(resp[4:], '%Y%m%d%H%M%S' )

                # keep only pdbs in the list that are new or newer, remove all others.
                if localdate == None or time.mktime(date) > time.mktime(localdate):
                    pdb_local_files[pdb] = date
                else:
                    requested_pdbs.remove(pdb)

            except error_perm:
                #550 error (permission error) results when a file does not exist, remove pdb from list
                requested_pdbs.remove(pdb)


        # ===============================================
        #3 download files that are in the list
        self.pdbTotal = len(requested_pdbs)
        self.processing_dates = False
        print 'Downloading %d PDB files' % self.pdbTotal

        for pdb in requested_pdbs:
            print "  [%d%%] - %s : %s" % (self.progress(), pdb, time.asctime(pdb_local_files[pdb])  )
            filename = 'pdb%s.ent.Z' % pdb
            local_filename = '%s/%s' % (ftp_update_settings.PDB_LOCAL_DIR,filename)

            #remove file if it exists already
            if os.path.exists(local_filename):
                os.remove(local_filename)

            self._incoming_file = open(local_filename,"w")
            ftp.retrbinary('RETR %s'%filename, self.processChunk, 1024)
            self._incoming_file.close()
            os.utime(local_filename, (time.mktime(pdb_local_files[pdb]),time.mktime(pdb_local_files[pdb])) )
            self.pdbCount += 1

        ftp.quit()
        print 'progress: %d%%' % self.progress()

    """
    Callback for ftp transfers.  This function will be called as chunks of data are received from the ftp server
    """
    def processChunk(self, data):
        self._incoming_file.write(data)


    """
    returns progress as a number between 0 and 100
    """
    def progress(self):
        if self.processing_dates:
            return 1
        elif self.pdbTotal == 0:
            return 100
        else:
            return round(self.pdbCount*1.0 / self.pdbTotal * 99)+1

        return self.pdbTotal

    """
    Returns the status as a string
    """
    def progressMessage(self):
        if self.processing_dates:
            return 'Processing list of PDBS'
        else:
            return '%d of %d PDB Files downloaded' % (self.pdbCount, self.pdbTotal)

    """
    Reset the task
    """
    def _reset(self):
        self.pdbCount = 0
        self.pdbTotal = 0


# code run if file executed from command line
if __name__ == '__main__':
    task = FTPUpdateTask('Command Line Update')

    args = {'3CK2A':1, '1FXKB':2, '1QE5A':3}

    task._work(args)




