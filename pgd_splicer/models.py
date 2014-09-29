""" ================================================================
# PDBSelect Selector Settings
================================================================ """
class PDBSelectSettings():
    RESOLUTION         = 1.75
    pdb_dir            = 'pdb'
    PDB_TMP_DIR        = 'tmp'
    DATA_VERSION       = '--'
pdb_select_settings = PDBSelectSettings()



""" ================================================================
# FTP Update Settings
================================================================ """
class FTPUpdateSettings():
    PDB_FTP_HOST   = 'ftp.ebi.ac.uk'
    PDB_REMOTE_DIR = '/pub/databases/rcsb/pdb/data/structures/all/pdb/'
    PDB_LOCAL_DIR  = './pdb/'
ftp_update_settings = FTPUpdateSettings()