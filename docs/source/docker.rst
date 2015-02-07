*********************
Using PGD with Docker
*********************

PGD ships with a Dockerfile to make development easier.
Consult the docker documentation for instructions how to use docker.
http://docs.docker.com/reference/


Using Fig
---------

Fig is a command line tool to automate using multiple docker containers.
Usually multiple long incantations of docker commands are necessary to get a
working development environment. Using fig, getting started takes just one
command:

::

   $ fig up

The following sections will not be necessary if you use fig.
Consult the fig documentation for details on how to modify the `fig.yml` file,
and other commands you can use with fig.
http://www.fig.sh/

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
and refer to using the docker command.

::

   $ docker run --name pgd_mysql \
    -e MYSQL_ROOT_PASSWORD=pgd \
    -e MYSQL_USER=pgd \
    -e MYSQL_PASSWORD=pgd \
    -e MYSQL_DATABASE=pgd \
    -d mysql

Running an Image and Linking it
-------------------------------

Once the MySQL container is running, we can run the PGD container we built and
link it with MySQL. Linking it means that the pgd container will be able to
transparently access it. We will also forward the container's port
8000 to the host's port 8000.

::

    $ docker run -d --name pgd -p 8000:8000 --link pgd_mysql:mysql osl_test/pgd

Mounting the PGD Code as a Volume
---------------------------------

Some developers may find the following to be convenient:

::

    $ docker run -d --name pgd \
    -p 8000:8000 \
    -v /path/to/code:/opt/pgd \
    --link pgd_mysql:mysql \
    osl_test/pgd
