*************************
Updating Protein Database
*************************

The PGD database is updated approximately four times per year. The following process should be followed to update the database.

One week before update: generate new selection file and download majority of updates.

    *From the pgd-staging Django site*::

        1 > python manage.py fetch --report=date-report.txt --selection=date-selection.txt

One day before update: reset staging database and load new selection file.

    *From the pgd-staging Django site*::

        1 > python manage.py shell


        1 > from pgd_core.models import Protein
        2 > for p in Protein.objects.all():
        3 >   p.delete()

        1 > ./pgd_splicer/ProcessPDBTask.py --pipein < date-selection.txt

On the update day: do the update!

    *From the pgd-staging Django site*::

        1 > ./pgd_splicer/dunbrack_selector.py --pipeout  > selection.txt
        2 > ./pgd_splicer/ftpupdate.py --pipein < selection.txt
        3 > ./pgd_splicer/ProcessPDBTask.py --pipein < selection.txt

That last command should be run within script so the output can be examined for particular failure modes once the command is complete.
Error messages and recommended actions:

    * "CRC check failed", "local variable 'i' referenced before assignment", "KeyError": collect all codes with these message, delete the files associated with these codes, retrieve them again from the server, and attempt to import the files again. Example command lines for codes 1B12 and 4AMW::

        1 > (cd ./pdb && for code in 1b12 4amw; do rm pdb$code.ent.gz; done)
        2 > egrep \(1B12\|4AMW\) selection.txt | ./pgd_splicer/ftpupdate.py --pipein
        3 > egrep \(1B12\|4AMW\) selection.txt | ./pgd_splicer/ProcessPDBTask.py --pipein

If the errors persist, document offending codes in bug reports as appropriate and copy the retrieved files aside for testing and comparison. * "Structure/DSSP mismatch", "No chains were parsed!": document offending codes in update post to PGD mailing list.

Cross-check database against selection file.

    *From the pgd-staging Django site*::

        1 > python manage.py crosscheck --selection=selection.txt

The staging site is now ready for customer preview.

Promote database from staging to production.

    *From the pgd-prod Django site*::

        1 > ./update-from-staging.sh

Rename the selection file for archival purposes. ::

    1 > mv selection.txt 201310-selection.txt

**On production, in settings.py, update the DATA_VERSION to the date used for the selection file.**
