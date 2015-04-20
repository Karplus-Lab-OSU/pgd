#!/bin/bash
EMPTY=0
HOWMANY=10
SHORTFILE=short-selection.txt

proteincount () {
    return 
}

# delete proteins from database
docker-compose run web bash -c 'echo "from pgd_core.models import Protein; Protein.objects.all().delete(); exit()" | python manage.py shell' >/dev/null

# the database should contain ${EMPTY} proteins
BEFORE=`docker-compose run web bash -c 'echo "from pgd_core.models import Protein; print Protein.objects.count(); exit()" | python manage.py shell' | tail -1 | awk '{ print $2 }' | tr -d '[[:space:]]'`
if [[ $BEFORE -ne $EMPTY ]]; then
    echo "FAIL: $((BEFORE)) proteins found, should be $((EMPTY))"
else
    echo "PASS: $((EMPTY)) proteins found"
fi

# generate selection file
docker-compose run web python ./pgd_splicer/dunbrack_selector.py --pipeout > selection.txt
sed ${HOWMANY}q selection.txt > ${SHORTFILE}

# retrieve files
docker-compose run web python ./pgd_splicer/ftpupdate.py --pipein < ${SHORTFILE}

# this command should list ${HOWMANY} files
FILES=`docker-compose run web ls -1 /opt/pgd/pdb | wc -l`
if [[ $FILES -ne $HOWMANY ]]; then
    echo "FAIL: $((FILES)) files found, should be $((HOWMANY))"
else
    echo "PASS: $((HOWMANY)) files found"
fi

# add proteins to database
docker-compose run web python ./pgd_splicer/ProcessPDBTask.py --pipein < ${SHORTFILE}

# the database should contain ${HOWMANY} proteins
AFTER=`docker-compose run web bash -c 'echo "from pgd_core.models import Protein; print Protein.objects.count(); exit()" | python manage.py shell' | tail -1 | awk '{ print $2 }' | tr -d '[[:space:]]'`
if [[ $AFTER -ne $HOWMANY ]]; then
    echo "FAIL: $((AFTER)) proteins found, should be $((HOWMANY))"
else
    echo "PASS: $((HOWMANY)) proteins found"
fi
