#!/usr/bin/env python
import os

from sets import Set
from ftplib import FTP
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

def processChunk(data):
    INCOMING_FILE.write(data)


pdb_remote_files = {}
pdb_ftpurl = 'ftp.ebi.ac.uk'
pdb_remote_dir = '/pub/databases/rcsb/pdb/data/structures/all/pdb/'

pdb_local_dir = '/var/www/pgd/splicer/pdb'
pdb_local_files = {}

timezone = 'pst'

#===============================================
#1 get list of local files matching the pdbs and the file modification date

#TODO For now simulate a list of pdbs being given to this function
#TODO the list will be created from the list of pdbs already in the directory
TEMP_local_files_all = os.listdir(pdb_local_dir + ".bak")
requested_pdbs = []
for filename in TEMP_local_files_all:
    requested_pdbs.append(filename[3:7])

for pdb in requested_pdbs:
    path = '%s/pdb%s.ent.Z' % (pdb_local_dir, pdb)
    if os.path.exists(path):
        date = time.gmtime(os.path.getmtime(path))
    else:
        date = None

    pdb_local_files[pdb] = date


print "LOCAL FILES:"
for entry in pdb_local_files:
    if pdb_local_files[entry] == None:
        print '  - %s : None' % (entry)
    else:
        print '  - %s : %s' % (entry, time.asctime(pdb_local_files[entry]))

print "------------------------------------------"

# ===============================================
#2 iterate through list and purge pdbs that are not new or newer
ftp = FTP(pdb_ftpurl)
#ftp.set_debuglevel(2)
ftp.login()
ftp.cwd(pdb_remote_dir)

for pdb in pdb_local_files:
    filename = 'pdb%s.ent.Z' % pdb

    localdate = pdb_local_files[pdb]

    resp = ftp.sendcmd('MDTM %s' % filename)
    #date = time.strptime('%s%s' % (resp[4:],time.tzname[0]), '%Y%m%d%H%M%S%Z' )
    date = time.strptime(resp[4:], '%Y%m%d%H%M%S' )

    # keep only pdbs that are new or newer, remove all others.
    if localdate == None or time.mktime(date) > time.mktime(localdate):
        pdb_local_files[pdb] = date
    else:
        requested_pdbs.remove(pdb)

# ===============================================
#3 download list of files
print 'Downloading %d PDB files' % len(requested_pdbs)

for pdb in requested_pdbs:
    print "  - %s : %s" % (pdb, time.asctime(pdb_local_files[pdb])  )
    filename = 'pdb%s.ent.Z' % pdb
    local_filename = 'tmp/%s' % filename

    #remove file if it exists already
    if os.path.exists(local_filename):
        os.remove(local_filename)

    INCOMING_FILE = open(local_filename,"w")
    ftp.retrbinary('RETR %s'%filename, processChunk, 1024)
    INCOMING_FILE.close()
    os.utime(local_filename, (time.mktime(pdb_local_files[pdb]),time.mktime(pdb_local_files[pdb])) )

ftp.quit()






