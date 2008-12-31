import dbsettings
from dbsettings.loading import set_setting_value

""" ================================================================
# PDBSelect Selector Settings
================================================================ """
class PDBSelectSettings(dbsettings.Group):
    RESOLUTION         = dbsettings.FloatValue('Resolution Limit', 'Import proteins with a resolution up to this limit')
    pdb_dir            = dbsettings.StringValue('PDB directory','Store pdbs in this directory')
    PDB_SELECT_URL     = dbsettings.StringValue('PDB Select URL', 'URL of directory where pdb select files are stored')
    PDB_TMP_DIR        = dbsettings.StringValue('PDB temp directory', 'Temporary directory for PDB files')
    PDB_SELECT_FILE_90 = dbsettings.StringValue('Threshold 90 Filename','Filename of PDB Select file containing proteins with a threshold of 90')
    PDB_SELECT_FILE_25 = dbsettings.StringValue('Threshold 25 Filename','Filename of PDB Select file containing proteins with a threshold of 25')
pdb_select_settings = PDBSelectSettings('Splicer')

# set defaults for PDBSelect settings
if not pdb_select_settings.RESOLUTION:
    set_setting_value('pgd_splicer.models', '', 'RESOLUTION', 1.75)
if not pdb_select_settings.pdb_dir:
    set_setting_value('pgd_splicer.models', '', 'pdb_dir', 'pdb')
if not pdb_select_settings.PDB_SELECT_URL:
    set_setting_value('pgd_splicer.models', '', 'PDB_SELECT_URL', 'http://bioinfo.tg.fh-giessen.de/pdbselect/')
if not pdb_select_settings.PDB_TMP_DIR:
    set_setting_value('pgd_splicer.models', '', 'PDB_TMP_DIR', 'tmp')
if not pdb_select_settings.PDB_SELECT_FILE_90:
    set_setting_value('pgd_splicer.models', '', 'PDB_SELECT_FILE_90', 'recent.pdb_select90')
if not pdb_select_settings.PDB_SELECT_FILE_25:
    set_setting_value('pgd_splicer.models', '', 'PDB_SELECT_FILE_25', 'recent.pdb_select25')


""" ================================================================
# FTP Update Settings
================================================================ """
class FTPUpdateSettings(dbsettings.Group):
    PDB_FTP_HOST   = dbsettings.StringValue('PDB FTP Host','Hostname of ftp server with PDB files')
    PDB_REMOTE_DIR = dbsettings.StringValue('PDB FTP Remote Directory','Directory on the ftp server where PDBs are stored')
    PDB_LOCAL_DIR  = dbsettings.StringValue('PDB Local Directory','Local Directory for storing PDB files')
ftp_update_settings = FTPUpdateSettings('Splicer')

# set defaults for FTP Update settings
if not ftp_update_settings.PDB_FTP_HOST:
    set_setting_value('pgd_splicer.models', '', 'PDB_FTP_HOST', 'ftp.ebi.ac.uk')
if not ftp_update_settings.PDB_REMOTE_DIR:
    set_setting_value('pgd_splicer.models', '', 'PDB_REMOTE_DIR', '/pub/databases/rcsb/pdb/data/structures/all/pdb/')
if not ftp_update_settings.PDB_LOCAL_DIR:
    set_setting_value('pgd_splicer.models', '', 'PDB_LOCAL_DIR', '/var/www/pgd/splicer/pdb')