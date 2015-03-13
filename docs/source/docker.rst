*********************
Using PGD with Docker
*********************

PGD ships with a Dockerfile to make development easier.
Consult the docker documentation for instructions on how to use docker.
http://docs.docker.com/reference/


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

