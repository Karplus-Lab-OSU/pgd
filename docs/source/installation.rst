************
Installation
************

===============================
Development Branch Installation
===============================

Dependencies::
    libmysqlclient-dev

In a virtualenv::
    pip install -r pgd/requirements.txt

==========================
Master Branch Installation
==========================
This is a manual installation guide for the Protein Geometry Database

^^^^^^^^^^^^^^^^^^^^
Install dependencies
^^^^^^^^^^^^^^^^^^^^

The following packages are required to run Protein Geometry Database. Either install these with your system's package manager or from `Pip <https://pip.pypa.io/en/latest/index.html>`_ (recommended). Note that system packages may not have the correct version.

    * `Python >= 2.5 but < 3.x <https://www.python.org/>`_ (Python 3.x is not supported due to backward-compatibility issues)
    * `setuptools >= 0.6.28 <https://pypi.python.org/pypi/setuptools>`_
    * `simplejson <https://pypi.python.org/pypi/simplejson>`_
    * `MySQL-python <https://pypi.python.org/pypi/MySQL-python>`_
    * `numpy <http://www.numpy.org/>`_
    * `biopython 1.57 <http://biopython.org/wiki/Main_Page>`_
    * `Python Django 1.3.x <https://docs.djangoproject.com/en/dev/intro/install/>`_
    * `Python Django-registration --0.7 <https://bitbucket.org/ubernostrum/django-registration/wiki/Home>`_
    * `Py2Cairo <http://cairographics.org/pycairo/>`_ See instructions below for installing the correct Py2Cairo
    * `Selenium <http://docs.seleniumhq.org/>`_ for running browser-based tests
    * liberation fonts - yum -y install git python-devel mysql-devel libffi-devel liberation-sans-fonts

The process for installing "py2cairo" was recently reported to be:

    #. Install libcairo2-dev via apt-get
    #. Install simplejson via pip to virtualenv
    #. Retrieve py2cairo from https://github.com/dieterv/py2cairo.git
    #. Run touch ChangeLog
    #. Run ./autogen.sh from that checkout
    #. Run ./waf configure --prefix=foo where foo is the absolute location of the virtualenv
    #. Run ./waf build
    #. Run ./waf install
    #. Open a python interpreter and confirm that the cairo module can be imported

^^^^^^^^^^^^
Get the Code
^^^^^^^^^^^^

    1. Make sure you have Git installed.
    2. Either download and unpack the "latest release":, or check it out from the repository::

        git clone https://github.com/osuosl/pgd

^^^^^^^^^^^^^
Configuration
^^^^^^^^^^^^^

    1. In the project root, you'll find a default-settings file called settings.py.dist. Copy it to settings.py::

        cp settings.py.dist settings.py

    2. If you want to use another database engine besides the default SQLite (not recommended for production), edit settings.py, and edit the following lines to reflect your wishes::

        1 DATABASE_ENGINE = ''   # <-- Change this to 'mysql'
        2 DATABASE_NAME = ''     # <-- Change this to a database name, or a file for SQLite
        3 DATABASE_USER = ''     # <-- Change this (not needed for SQLite)
        4 DATABASE_PASSWORD = '' # <-- Change this (not needed for SQLite)
        5 DATABASE_HOST = ''     # <-- Change this (not needed if database is localhost)
        6 DATABASE_PORT = ''     # <-- Change this (not needed if database is localhost)

    3. Initialize Database::

        ./manage.py syncdb

    4. Everything should be all set up! Run the development server with::

        ./manage.py runserver

--------------
Importing Data
--------------

See the `importing data
<https://code.osuosl.org/projects/pgd/wiki/Designsplicercli>`_ section for instructions on how to import data.

------------------------------------------------
Additional configuration for production servers:
------------------------------------------------

Deploying a production server requires additional setup steps.

    1. Change your SECRET_KEY to unique (and hopefully unguessable) strings in your settings.py.
    2. Ensure the server has the ability to send emails or you have access to an SMTP server. Set EMAIL_HOST, EMAIL_PORT, and DEFAULT_FROM_EMAIL in settings.py. For more complicated outgoing mail setups, please refer to the django email documentation.
    3. Follow the django guide to deploy with apache. Here is an example mod_wsgi file::

        1  import os
        2  import sys
        3
        4  path = '/var/lib/django/pgd'
        5  if path not in sys.path:
        6      sys.path.append(path)
        7
        8  os.environ['DJANGO_SETTINGS_MODULE'] = 'settings'
        9
        10 import django.core.handlers.wsgi
        11 application = django.core.handlers.wsgi.WSGIHandler()

    4. If New Relic support is required, modify the WSGI file according to the New Relic documentation. On zeus, the New Relic configuration file is /etc/newrelic.ini.
