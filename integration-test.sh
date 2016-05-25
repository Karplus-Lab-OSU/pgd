#!/bin/bash

# takes an optional argument: how many proteins to retrieve and import

EMPTY=0
HOWMANY=${1:-10}

SELFILE=$(mktemp)
SHORTFILE=$(mktemp)

# NB: these tests only run on docker!

# check if we are already running on docker
if [ ! -f /.dockerenv ]; then
    DOCKER='docker-compose run --rm web'
else
    DOCKER=
fi

# delete proteins from database
$DOCKER python manage.py coredb --clear

# the database should contain ${EMPTY} proteins
TESTEMPTY=0
BEFORE=`$DOCKER python manage.py coredb --count | grep proteins | awk '{ print $1 }'`
if [[ $BEFORE -ne $EMPTY ]]; then
    echo "FAIL: $((BEFORE)) proteins found, should be $((EMPTY))"
else
    echo "PASS: $((EMPTY)) proteins found"
    TESTEMPTY=1
fi

# generate selection file
$DOCKER python ./pgd_splicer/dunbrack_selector.py --pipeout > ${SELFILE}
sed ${HOWMANY}q ${SELFILE} > ${SHORTFILE}

# retrieve files
$DOCKER find /opt/pgd/pdb -type f -exec rm {} \;
$DOCKER python ./pgd_splicer/ftpupdate.py --pipein < ${SHORTFILE}

# this command should list ${HOWMANY} files
TESTFILES=0
FILES=`$DOCKER find /opt/pgd/pdb -type f -print | wc -l`
if [[ $FILES -ne $HOWMANY ]]; then
    echo "FAIL: $((FILES)) files found, should be $((HOWMANY))"
else
    echo "PASS: $((HOWMANY)) files found"
    TESTFILES=1
fi

# add proteins to database
$DOCKER python ./pgd_splicer/ProcessPDBTask.py --pipein < ${SHORTFILE}

# the database should contain ${HOWMANY} proteins
TESTFULL=0
AFTER=`$DOCKER python manage.py coredb --count | grep proteins | awk '{ print $1 }'`
if [[ $AFTER -ne $HOWMANY ]]; then
    echo "FAIL: $((AFTER)) proteins found, should be $((HOWMANY))"
else
    echo "PASS: $((HOWMANY)) proteins found"
    TESTFULL=1
fi

rm -rf ${SELFILE} ${SHORTFILE}

# final score
if [[ ${TESTEMPTY} -eq 1 && ${TESTFILES} -eq 1 && ${TESTFULL} -eq 1 ]]; then
    exit 0
else
    exit 1
fi
