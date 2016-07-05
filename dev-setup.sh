#!/bin/bash

# clean up current situation
docker-compose stop
docker-compose rm -v --force --all
find . \( -name "*.pyc" -o -name "*.pyo" \) -print0 | xargs -0 rm -rf

# rebuild containers
docker-compose build

# start data containers
docker-compose up -d mysql

# collectstatic
docker-compose run --rm web python manage.py collectstatic --noinput

# create database
# due to race conditions, this may require multiple attempts
retcode=1
count=0
syncdbmax=10
while [ $retcode -ne 0 -o $count -eq $syncdbmax ]; do
    docker-compose run --rm web python manage.py syncdb --noinput >/dev/null 2>/dev/null
    count=`expr $count + 1`
    retcode=$?
done
if [ $count -eq $syncdbmax ]; then
    echo "syncdb failed:"
    docker-compose run --rm web python manage.py syncdb --noinput
fi
