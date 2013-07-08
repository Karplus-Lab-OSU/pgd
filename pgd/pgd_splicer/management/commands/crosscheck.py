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
        self.infile = []

        # Dict of proteins in database with settings that differ from
        # those in the selection file.
        # Key = protein code
        # Value = Tuple of two elements: file and database settings
        self.wrong = {}

        # Dict of proteins in database with missing chains.
        # Key = protein code
        # Value = List of missing chains
        self.chains = {}

        # for every protein in the selection file:
        #   is the protein in the database?
        #   -- no, add to notfound list and continue
        #   are the settings correct?
        #   -- no, add to wrongsettings dict and continue
        #   for every chain in the protein
        #     is the chain in the database?
        #     -- no, add to nochains dict and continue

        # Iterate through selection file.
        with open(self.selection, 'r') as f:
            line = f.readline()
            try:
                code, selchains, threshold, resolution, rfactor, rfree = line
            except ValueError:
                # Version line only has two entries, we can skip it
                continue

            # Append code to list found in file
            self.infile.append(code)

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

        # Proteins found in one place but not in the other.
        self.notfile = list(self.indb - self.infile)
        self.notdb = list(self.infile - self.indb)

        # Generate report.
        if self.notfile is []:
            sys.stdout.write('All proteins in the database are in the file.')
        else:
            sys.stdout.write('Proteins found in the database ',
                             'but not in the file: %d\n' % len(self.notfile))
            sys.stdout.write(', '.join(sorted(self.notfile)), '\n')

        if self.notdb is []:
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
