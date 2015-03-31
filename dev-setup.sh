#!/bin/bash

# clean up current situation
fig stop
fig rm --force

# rebuild containers
fig build

# start data containers
fig up -d mysql pdb

# create database
# due to race conditions, this may require multiple attempts
retcode=1
while [ $retcode -ne 0 ]; do
    fig run web python manage.py syncdb --noinput >/dev/null 2>/dev/null
    retcode=$?
done
