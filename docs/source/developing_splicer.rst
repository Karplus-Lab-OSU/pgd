******************
Developing Splicer
******************

These are some special instructions for developing splicer components.

---------------
Django Settings
---------------

Django and its ORM can be used outside of its webserver. The only requirements are
    * the directory containing **settings.py** is on the sys.path
    * environment variable **DJANGO_SETTINGS_MODULE** is set to the **settings**

Splicer components should all do this automatically using python code to add the correct directory to **sys.path**. This works as long as the components are run from the directory containing **settings.py** ::

    import os, sys
    #python magic to add the current directory to the pythonpath
    sys.path.append(os.getcwd())

    # ----------------------------------------------------------
    # Setup django environment 
    # ----------------------------------------------------------
    if not os.environ.has_key('DJANGO_SETTINGS_MODULE'):
    os.environ['DJANGO_SETTINGS_MODULE'] - 'settings
    # ----------------------------------------------------------
    # Done setting up django environment
    # ----------------------------------------------------------

----------------------------------------
Running Components from the command line
----------------------------------------

Splicer Subtasks all contain **main** code that starts the tasks using arguments passed in. All should give you a list of required arguments if you do not pass them in. The main intention of this is to allow easier debugging on tasks.

For example **ProcessPDBTask** requires a pdb code followed by resolution, threshold, rfactor, and rfree. For testing purposing the properties do not need to be real values. ::

    ~/pgd/pgd $ ./pgd_splicer/ProcessPDBTask.py 12SB 1 2 3 4
