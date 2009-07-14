import dbsettings


""" ================================================================
# PDBSelect Selector Settings
================================================================ """
class PDBSelectSettings(dbsettings.Group):
    RESOLUTION         = dbsettings.FloatValue('Resolution Limit', 'Import proteins with a resolution up to this limit', default=1.75)
    pdb_dir            = dbsettings.StringValue('PDB directory','Store pdbs in this directory', default='pdb')
    PDB_SELECT_URL     = dbsettings.StringValue('PDB Select URL', 'URL of directory where pdb select files are stored',default='http://bioinfo.tg.fh-giessen.de/pdbselect/')
    PDB_TMP_DIR        = dbsettings.StringValue('PDB temp directory', 'Temporary directory for PDB files','tmp')
    PDB_SELECT_FILE_90 = dbsettings.StringValue('Threshold 90 Filename','Filename of PDB Select file containing proteins with a threshold of 90',default='recent.pdb_select90')
    PDB_SELECT_FILE_25 = dbsettings.StringValue('Threshold 25 Filename','Filename of PDB Select file containing proteins with a threshold of 25',default='recent.pdb_select25')
pdb_select_settings = PDBSelectSettings('Splicer')



""" ================================================================
# FTP Update Settings
================================================================ """
class FTPUpdateSettings(dbsettings.Group):
    PDB_FTP_HOST   = dbsettings.StringValue('PDB FTP Host','Hostname of ftp server with PDB files', default='ftp.ebi.ac.uk')
    PDB_REMOTE_DIR = dbsettings.StringValue('PDB FTP Remote Directory','Directory on the ftp server where PDBs are stored', default='/pub/databases/rcsb/pdb/data/structures/all/pdb/')
    PDB_LOCAL_DIR  = dbsettings.StringValue('PDB Local Directory','Local Directory for storing PDB files', default='./pdb/')
ftp_update_settings = FTPUpdateSettings('Splicer')