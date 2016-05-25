from django.core.management.base import BaseCommand, CommandError
from optparse import make_option
from pgd_core.models import Protein, Chain, Residue
# from django.conf import settings
# from pgd_splicer.tools import localfile, remotefile
# import urllib
# import re
# import gzip
# from cStringIO import StringIO
# import os
# import time
# from ftplib import FTP, error_perm
# from datetime import datetime


class Command(BaseCommand):
    option_list = BaseCommand.option_list + (
        make_option('--count',
                    action='store_true',
                    dest='count',
                    default=False,
                    help='Count proteins, chains, and residues in core ' +
                    'database.'),
        make_option('--clear',
                    action='store_true',
                    dest='clear',
                    default=False,
                    help='Clear core database'),
    )
    help = 'Performs basic functions on core database.'

    methods = {}

    def countdb(self):
        # Count proteins, chains, and residues in core database.
        self.stdout.write('{} proteins'.format(Protein.objects.all().count()))
        self.stdout.write('{} chains'.format(Chain.objects.all().count()))
        self.stdout.write('{} residues'.format(Residue.objects.all().count()))
    methods['count'] = countdb

    def cleardb(self):
        # Clear core database of all proteins
        for p in Protein.objects.all():
            p.delete()
    methods['clear'] = cleardb

    def handle(self, *args, **options):
        # Only one value is allowed in options
        possible = [option for option in options if options[option] is True]
        if len(possible) != 0:
            CommandError('Exactly one option must be enabled')
        option = possible[0]

        # If value in options is not in methods, complain.
        if option not in self.methods:
            CommandError('Option {} not in methods'.format(option))

        # Execute method.
        self.methods[option](self)
