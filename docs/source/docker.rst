*********************
Using PGD with Docker
*********************

PGD ships with a Dockerfile to make development easier.
Consult the docker documentation for instructions on how to use docker.
http://docs.docker.com/reference/


Quick Start: Demonstrating the PGD with Fig and Docker
------------------------------------------------------

This is for folks who are already familiar with Fig and Docker.  If
you are new to either of these tools, please skip ahead to the next
section of the documentation.

The first step is building the containers.


::

   $ fig build

Now create the necessary database tables.  For this version of Django,
syncdb is still required.  The PGD does not use the admin site at this
time, so there's no need to create an account.

::

   $ fig run web python manage.py syncdb --noinput

This command may fail the first time with a lack of connection due to
fig's not-yet-mature orchestration functionality.  Simply run it again
and it should succeed.

The next step is to create a selection file.  It's possible to use a
subset of an existing selection file (say the first hundred lines) but
if you need to generate a new one, use this command:

::

   $ fig run web python ./pgd_splicer/dunbrack_selector.py --pipeout > selection.txt
   $ sed 100q selection.txt > top-100-selection.txt

The selected proteins must be retrieved from the worldwide PDB
collection.  This command may take some time!

::

   $ fig run web python ./pgd_splicer/ftpupdate.py --pipein < top-100-selection.txt

To list the proteins that were successfully downloaded, run this command:

::

   $ fig run web ls /opt/pgd/pdb

You should see 100 files with names like `pdb1ae1.ent.gz`.

Finally, all the retrieved proteins must be imported into the
database.  This command will definitely take some time: a full update
currently consists of over twenty-six thousand proteins, and can take
upwards of eight hours to process.

::

   $ fig run web python ./pgd_splicer/ProcessPDBTask.py --pipein < top-100-selection.txt

To confirm the number of proteins in the database, use the Django shell:

::

   $ fig run web python manage.py shell
   Python 2.7.5 (default, Jun 17 2014, 18:11:42) 
   [GCC 4.8.2 20140120 (Red Hat 4.8.2-16)] on linux2
   Type "help", "copyright", "credits" or "license" for more information.
   (InteractiveConsole)
   >>> from pgd_core.models import Protein
   >>> Protein.objects.count()
   100
   >>> 

Looking good!  Now it's time to actually start the web server.

::

   $ fig up

This will generate a screen or two of output from the different
containers.  Once that output stabilizes, open a web browser to
`http://localhost:8000` (or a different host, depending on where
you're running Docker) and you should see the PGD!  Select 'Search',
remove the default search constraints on omega from the search page,
and select 'Submit', and you should see a Ramachandran plot with
results.  Success!

Using Fig
---------

Fig is a command line tool to automate using multiple docker
containers.  Usually multiple long incantations of docker commands are
necessary to get a working development environment. Using fig, a
simple instance of the PGD without any content can be started from
scratch with a simple command:

::

   $ fig up

The application will be available on ``http://localhost:8000``.

.. note::
	Fig is a developer tool, and this early version is prone to certain kinds
	of race conditions. It is possible for the web container to come up before
	the database container, and if the web container can't find the database it
	will fail.

.. note::
	Fig has been renamed docker-compose. Our workflow should be updated to
	reflect this.

	- The fig command will be renamed docker-compose
	- The fig.yml file will be renamed to .docker-compose.yml
	- The PyPI package will be renamed to docker-compose
	- These docs will need to be updated.

Similarly, to run all the tests in the PGD code base, the following
command can be very useful:

::

   $ fig run web python manage.py test

Consult the fig documentation for details on how to modify the `fig.yml` file,
and other commands you can use with fig.
http://www.fig.sh/
   
The following sections will not be necessary if you use fig.

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

