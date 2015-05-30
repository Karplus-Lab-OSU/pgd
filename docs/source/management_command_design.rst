*****************************************
Management Commands and Possible Redisign
*****************************************

2 Databases:
    * Staging (aka Silver)
    * Master (aka Gold)

Instead of two databases ::

    class Protein(models.Model):
        """ 
        Same as before
        """ 
        ...

    class GoldProtein(Protein):
        """ 
        Nothing actually goes here
        """ 

Manage Cammonds:

    * Import (Modifies Staging, Reads from Master)
        * Fetches pdb files (like fetch does currently)
            * --fetch-only as an option
        * Stores the selection in the Audit table
        * Dumps Proteins from staging
        * "ProcessPDBTask"
        * Generates a diff (Total, New, Removed) data and stores it in the Audit table

    * Promote (Reads Staging, Modifies Master)
        * Dumps Data
        * Copies Staging into Master
        * Updates Master Audit table

