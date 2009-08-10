import dbsettings


""" ================================================================
# PDBSelect Selector Settings
================================================================ """
class PDBSelectSettings(dbsettings.Group):
    RESOLUTION         = dbsettings.FloatValue('Resolution Limit', 'Import proteins with a resolution up to this limit', default=1.75)
    pdb_dir            = dbsettings.StringValue('PDB directory','Store pdbs in this directory', default='pdb')
    PDB_TMP_DIR        = dbsettings.StringValue('PDB temp directory', 'Temporary directory for PDB files','tmp')
    DATA_VERSION       = dbsettings.StringValue('Version of imported data', 'This is something to identify the version of file that splicer imported data from.  This may be an actual version number, or the date the select file was generated', default='--')
pdb_select_settings = PDBSelectSettings('Splicer')



""" ================================================================
# FTP Update Settings
================================================================ """
class FTPUpdateSettings(dbsettings.Group):
    PDB_FTP_HOST   = dbsettings.StringValue('PDB FTP Host','Hostname of ftp server with PDB files', default='ftp.ebi.ac.uk')
    PDB_REMOTE_DIR = dbsettings.StringValue('PDB FTP Remote Directory','Directory on the ftp server where PDBs are stored', default='/pub/databases/rcsb/pdb/data/structures/all/pdb/')
    PDB_LOCAL_DIR  = dbsettings.StringValue('PDB Local Directory','Local Directory for storing PDB files', default='./pdb/')
ftp_update_settings = FTPUpdateSettings('Splicer')