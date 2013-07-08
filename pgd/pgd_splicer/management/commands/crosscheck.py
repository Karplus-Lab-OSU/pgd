from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from pgd_core.models import Protein, Chain
import sys
import os


class Command(BaseCommand):
    option_list = BaseCommand.option_list + (
        make_option('--selection',
                    default=False,
                    help='write selections to FILE'),
    )
    help = 'Cross-checks database against selection file.'

    def handle(self, *args, **options):
        self.selection = options['selection']
        if not self.selection:
            raise CommandError('Selection file required!')
        if not os.path.exists(self.selection):
            raise CommandError('Selection file does not exist')

        # List of protein codes in database.
        self.indb = set(Protein.objects.all().values_list('code', flat=True))

        # List of proteins found in file.
        self.infile = set([])

        # Dict of proteins in database with settings that differ from
        # those in the selection file.
        # Key = protein code
        # Value = Tuple of two elements: file and database settings
        self.wrong = {}

        # Dict of proteins in database with missing chains.
        # Key = protein code
        # Value = List of missing chains
        self.chains = {}

        # Iterate through selection file.
        f = open(self.selection, 'r')
        for line in f:
            try:
                code, selchains, threshold, resolution, rfactor, rfree = line
            except ValueError:
                # Version line only has two entries, we can skip it
                continue

            # Append code to list found in file
            self.infile.update(code)

            # If protein found in database, perform remaining checks.
            if code in self.indb:
                protein = Protein.objects.get(code=code)
                fileset = (code, threshold, resolution, rfactor, rfree)
                dbset = (protein.code, protein.threshold, protein.resolution,
                         protein.rfactor, protein.rfree)
                if fileset != dbset:
                    self.wrong[code] = (fileset, dbset)

                # Check chains.
                badchains = []
                for selchain in list(selchains):
                    try:
                        chainid = '%s%s' % (code, selchain)
                        chain = protein.chains.get(id=chainid)
                    except Chain.DoesNotExist:
                        badchains.append(selchain)
                if badchains is not []:
                    self.chains[code] = badchains
        f.close()

        # Proteins found in one place but not in the other.
        self.notfile = self.indb.difference(self.infile)
        self.notdb = self.infile.difference(self.indb)

        # Generate report.
        if len(self.notfile) == 0:
            sys.stdout.write('All proteins in the database are in the file.')
        else:
            sys.stdout.write('Proteins found in the database ',
                             'but not in the file: %d\n' % len(self.notfile))
            sys.stdout.write(', '.join(sorted(self.notfile)), '\n')

        if len(self.notdb) == 0:
            sys.stdout.write('All proteins in the file are in the database.')
        else:
            sys.stdout.write('Proteins found in the file ',
                             'but not in the database: %d\n' % len(self.notdb))
            sys.stdout.write(', '.join(sorted(self.notdb)), '\n')

        if self.wrong is {}:
            sys.stdout.write('All proteins in the database ',
                             'have settings which match the selection file.')
        else:
            sys.stdout.write('Proteins found in the database ',
                             'with incorrect settings: %d\n' % len(self.wrong))
            for key in sorted(self.wrong.iterkeys()):
                sys.stdout.write('  %s: %s\n' % (key, self.wrong[key]))

        if self.chains is {}:
            sys.stdout.write('All proteins in the database ',
                             'have all chains listed in the selection file.')
        else:
            sys.stdout.write('Proteins found in the database ',
                             'missing chains: %d\n' % len(self.chains))
            for key in sorted(self.chains.iterkeys()):
                sys.stdout.write('  %s: %s\n' % (key, self.chains[key]))
