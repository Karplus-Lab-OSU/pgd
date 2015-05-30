************
Installation
************

This is a manual installation guide for the Protein Geometry Database (PGD).

^^^^^^^^^^^^^^^^^^^^
Install dependencies
^^^^^^^^^^^^^^^^^^^^

Certain operating-system packages are required to run the PGD.  As an
example, here are the additional packages required for a Centos 7
server with EPEL:

* bzip2 
* cairo-devel 
* gcc 
* gcc-c++ 
* libffi 
* libffi-devel 
* mysql 
* mysql-devel 
* nodejs 
* npm 
* python-devel 
* python-setuptools 
* tar

Also required is ``dsspcmbi``, compiled from the DSSP_ software.

.. _DSSP: http://swift.cmbi.ru.nl/gv/dssp/

Python 2.7 or greater is required, but Python 3.x is not currently supported.

^^^^^^^^^^^^
Get the Code
^^^^^^^^^^^^

#. Make sure you have Git installed.

#. Check it out from the repository::

     git clone https://github.com/osuosl/pgd

^^^^^^^^^^^^^
Configuration
^^^^^^^^^^^^^

The easiest way to spin up an instance of the PGD for testing purposes
is with :doc:`Docker <docker>`.  If Docker is unavailable, an instance
of the PGD can be spun up locally following these instructions.

#. Construct a ``settings.ini`` file in the top level of the
   application.  This file must contain values for ``SECRET_KEY``,
   with other values optional.  Here is an example::

     [settings]
     MYSQL_ENV_MYSQL_DATABASE=pgd
     MYSQL_ENV_MYSQL_USER=pgd
     MYSQL_ENV_MYSQL_PASSWORD=kjwb_if4hgkujpb3*7(_8
     MYSQL_PORT_3306_TCP_ADDR=mysql.example.org
     MYSQL_PORT_3306_TCP_PORT=3306
     GOOGLE_ID=UA-8675309-1
     SECRET_KEY=2nWoBbgLb1bVbOzM0PaW/q0jScKKcP5j2nWoBbgLb1bVbOzM0P
     MEDIA_ROOT=/opt/django/pgd/media
     STATIC_ROOT=/opt/django/pgd/static
     EMAIL_HOST=smtp.example.org
     DEFAULT_FROM_EMAIL=registration@pgd.example.org
     SERVER_EMAIL=pgd@pgd.example.org

#. Initialize the database::

     $ python manage.py syncdb

#. Collect the static files::
	  
     $ python manage.py collectstatic

#. Now the server can be run::

     $ python manage.py runserver

See the :doc:`importing data <importing_data>` section for instructions on how to import data.

More information about deploying with WSGI can be found here_.

.. _here: https://docs.djangoproject.com/en/1.6/howto/deployment/wsgi/
