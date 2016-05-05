from django.conf import settings

import os


# localfile returns the location of the PDB file corresponding to the given code.
def localfile(code):
    lc = code.lower()
    subdir = lc[1:3]
    name = 'pdb{}.ent.gz'.format(lc)
    return os.path.join(settings.PDB_LOCAL_DIR, subdir, name)


# remotefile does the same for remote files
def remotefile(code):
    lc = code.lower()
    subdir = lc[1:3]
    name = 'pdb{}.ent.gz'.format(lc)
    return os.path.join(settings.PDB_REMOTE_DIR, subdir, name)
