#!/bin/bash

# clean up current situation
docker-compose stop
docker-compose rm --force

# rebuild containers
docker-compose build

# start data containers
docker-compose up -d mysql pdb

# collectstatic
docker-compose run web python manage.py collectstatic --noinput

# create database
# due to race conditions, this may require multiple attempts
retcode=1
while [ $retcode -ne 0 ]; do
    docker-compose run web python manage.py syncdb --noinput >/dev/null 2>/dev/null
    retcode=$?
done
