*********************
Using PGD with Docker
*********************

PGD ships with a Dockerfile to make development easier.
Consult the docker documentation for instructions on how to use docker.
http://docs.docker.com/reference/


Quick Start: Demonstrating the PGD with Docker Compose and Docker
-----------------------------------------------------------------

This is for folks who are already familiar with Docker Compose and
Docker.  If you are new to either of these tools, please skip ahead to
the next section of the documentation.

There are three steps to this process: building the containers,
populating them with content, and starting the web server.

Build the containers
====================

.. note::

   The repository contains a script named `dev-setup.sh` which builds
   the containers following the same instructions found in this
   script.  Use at your own risk as the script may not be updated as
   often as the documentation.  When in doubt, trust the docs!
   
This is pretty straightforward.


::

   $ docker-compose build

The database containers (for MySQL and PDB files) need to be brought
up next.

::

   $ docker-compose up -d mysql pdb

   
Now create the necessary database tables.  For this version of Django,
syncdb is still required.  The PGD does not use the admin site at this
time, so there's no need to create an account.

::

   $ docker-compose run web python manage.py syncdb --noinput

This command may fail the first time with a lack of connection due to
docker-compose's not-yet-mature orchestration functionality.  Simply
run it again and it should succeed.

Install the content
===================

.. note::

   The repository contains a script named `integration_test.sh` which
   tests the scripts mentioned in this section.  Use at your own risk
   as the script may not be updated as often as the documentation.
   Again, when in doubt, trust the docs!

The next step is to create a selection file.  It's possible to use a
subset of an existing selection file (say the first hundred lines) but
if you need to generate a new one, use this command:

::

   $ docker-compose run web python ./pgd_splicer/dunbrack_selector.py --pipeout > selection.txt
   $ sed 100q selection.txt > top-100-selection.txt

The selected proteins must be retrieved from the worldwide PDB
collection.  This command may take some time!

::

   $ docker-compose run web python ./pgd_splicer/ftpupdate.py --pipein < top-100-selection.txt

To list the proteins that were successfully downloaded, run this command:

::

   $ docker-compose run web ls /opt/pgd/pdb

You should see 100 files with names like `pdb1ae1.ent.gz`.

Finally, all the retrieved proteins must be imported into the
database.  This command will definitely take some time: a full update
currently consists of over twenty-six thousand proteins, and can take
upwards of eight hours to process.

::

   $ docker-compose run web python ./pgd_splicer/ProcessPDBTask.py --pipein < top-100-selection.txt

To confirm the number of proteins in the database, use the Django shell:

::

   $ docker-compose run web python manage.py shell
   Python 2.7.5 (default, Jun 17 2014, 18:11:42) 
   [GCC 4.8.2 20140120 (Red Hat 4.8.2-16)] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
   (InteractiveConsole)
   >>> from pgd_core.models import Protein
   >>> Protein.objects.count()
   100
   >>> 

Start it up!
============   
   
Looking good!  Now it's time to actually start the web server.

::

   $ docker-compose up

This will generate a screen or two of output from the different
containers.  Once that output stabilizes, open a web browser to
`http://localhost:8000` (or a different host, depending on where
you're running Docker) and you should see the PGD!  Select 'Search',
remove the default search constraints on omega from the search page,
and select 'Submit', and you should see a Ramachandran plot with
results.  Success!

Using Docker Compose
--------------------

Docker Compose is a command line tool to automate using multiple
docker containers.  Usually multiple long incantations of docker
commands are necessary to get a working development environment. Using
docker-compose, a simple instance of the PGD without any content can
be started from scratch with a simple command:

::

   $ docker-compose up

The application will be available on ``http://localhost:8000``.

.. note::
	Docker-Compose is a developer tool, and this early version is
	prone to certain kinds of race conditions. It is possible for
	the web container to come up before the database container,
	and if the web container can't find the database it will fail.

Similarly, to run all the tests in the PGD code base, the following
command can be very useful:

::

   $ docker-compose run web python manage.py test

Consult `the docker-compose documentation
<http://docs.docker.com/compose/>`_ for details on how to modify the
`docker-compose.yml` file, and other commands you can use with
docker-compose.
   
The following sections will not be necessary if you use docker-compose.

Building an Image
-----------------

To build an image with PGD installed, run this command:

::

   $ docker build -t osl_test/pgd .

The `-t` option specifies the tag for the image. We use `osl_test` here for
testing.

Running a MySQL Container
-------------------------

PGD relies on a MySQL database. We use the default `mysql` image. Docker will
fetch the `mysql` image automatically.
The `-e` option passes an environment variable to the image. In this example we
set a series of necessary environment variables to a simple default.
The `--name` option gives this new container a name so it is easier to remember
and reference when using the docker command.

::

   $ docker run --name pgd_mysql \
    -e MYSQL_ROOT_PASSWORD=pgd_root_password \
    -e MYSQL_USER=pgd_user \
    -e MYSQL_PASSWORD=pgd_user_password \
    -e MYSQL_DATABASE=pgd_db \
    -d mysql

Running an Image and Linking it
-------------------------------

Once the MySQL container is running, we can run the PGD container we built and
link it with MySQL. Linking it means that the pgd container will be able to
transparently access it. We will also forward the container's port
8000 to the host's port 8000.

::

    $ docker run -d --name pgd -p 8000:8000 --link pgd_mysql:mysql osl_test/pgd

This should result in an instance of the PGD running on localhost at port 8000.
       
Mounting the PGD Code as a Volume
---------------------------------

Some developers may find the following to be convenient:

::

    $ docker run -d --name pgd \
    -p 8000:8000 \
    -v /path/to/code:/opt/pgd \
    --link pgd_mysql:mysql \
    osl_test/pgd

Be warned: this may clash with the Dockerfile's treatment of
`settings.py` depending on whether one already exists in the checkout.

