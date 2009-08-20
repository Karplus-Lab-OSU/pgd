""" *******************************************************************************
    Prints out the query results in a HTML format or in a dump file
    in an HTML format the function will paginate the data
    numResidues = total # of residues in query
    Result_Set = set of results to be printed
    Per_Page = number of results to show per page
    toFile = file to dump results to
    **************************************************************************  """
from __future__ import with_statement
from threading import Thread, Lock


import math
from pgd_search.models import *
from pgd_constants import AA_CHOICES
from pgd_search.models import searchSettings
from django.core.paginator import Paginator


# A list of values that should not be printed out
FIELDS = ['aa','a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','ss','phi', 'psi', 'ome', 'chi', 'bm', 'bs', 'bg', 'h_bond_energy', 'zeta']
FIELD_LABEL_REPLACEMENTS = {
    'h_bond_energy':'H Bond', 
    'aa':'AA',
    'L1':u'C(-1)N',
    'L2':u'N-CA',
    'L3':u'CA-CB',
    'L4':u'CA-C',
    'L5':'C-O',
    'a1':u'C(-1)-N-CA',
    'a2':u'N-CA-CB',
    'a3':u'N-CA-C',
    'a4':u'CB-CA-C',
    'a5':u'CA-C-O',
    'a6':u'CA-C-N(+1)',
    'a7':u'O-C-N(+1)'
    }
FIELD_VALUE_REPLACEMENTS = {'aa':AA_CHOICES}


class BufferThread(Thread):
    """
    Adds processed segments to a Dump's buffer.  This class
    is used so that the buffering happens at the same time as the file
    being read.

    Page sizes are defined by the Dump class.  Each page contains a set
    of segments which are processed into text lines and added to the buffer
    """

    def __init__(self, parent):
        print 'bufferthread created'
        self.parent = parent
        Thread.__init__(self)

    def run(self):
        self.parent

        lines = []

        page_num = self.parent.current_page
        self.parent.current_page += 1
        for segment in self.parent.pages.page(page_num).object_list:
            self.parent.count += 1
            first = True
            for offset, string in self.parent.iValues:

                parts = [
                    str(self.parent.count) if first else ' ',
                    segment.protein_id,
                    string,
                    segment.__dict__['r%i_oldID' % (self.parent.iIndex+offset) ],
                    segment.chainID,
                ]

                #field values
                for field in FIELDS:
                    # replace field with display value if needed
                    if field in FIELD_VALUE_REPLACEMENTS:
                        code = segment.__dict__['r%i_%s' % (iIndex+offset, field)]
                        if code:
                            for k,v in FIELD_VALUE_REPLACEMENTS[field]:
                                if k == code:
                                    parts.append(str(v))
                    # just write value
                    else:
                        parts.append(str(segment.__dict__['r%i_%s' % (iIndex+offset, field)]))

                s = '\t'.join(parts)
                string = '%s\n' % s

                with self.parent.buffer_lock:
                    self.parent.buffer.append(string)

        # update parents buffer with new lines then release this thread so 
        # parent can create a new thread if needed
        with self.parent.buffer_lock:
            if self.parent.current_page > self.parent.page_max:
                self.parent.nEOF = False

            # dereference self so that another thread can run
            self.parent.buffer_thread = None


class Dump():
    """
    Class that encapsulates a dump of a queryset.  This class turns the 
    results of the query set into an iterable returning sections of text
    that make up the dump file.
    """

    # maximum size for the buffer
    buffer_size = 20000

    # number of rows to increment the buffer by
    # keep in mind that each row will generate multiple lines
    # depending on how big the segment_length is
    buffer_increment = 500


    def __init__(self, search):

        self.buffer = []
        self.buffer_thread = None
        self.buffer_lock = Lock()

        self.pages = Paginator(search.querySet(), self.buffer_increment)
        self.page_max = self.pages.page_range[-1]
        self.current_page = 1
        self.count = 0
        self.nEOF = True
        self.create_header()

        #calculate list of iValues
        self.iValues = [
            (
                i, #index int
                ('(i%+i)' % i) if i else '(i)', # string representation
            ) for i in range(
                0 - (search.segmentLength-1)/2, #start
                int(math.ceil((search.segmentLength-1) / 2.0))+1, #stop
            )
        ] 

        #calculate iIndex
        self.iIndex = int(math.ceil(searchSettings.segmentSize/2.0)-1)


    def create_header(self):
        """
        Creates the header for the dump.  This must happen before any calls to
        next otherwise the update threads will begin filling the buffer with
        output from the rows in the queryset.
        """

        parts = ["Match\tCode\tResidue\tID\tChain ID"]

        # Field names
        for field in FIELDS:
            if field in FIELD_LABEL_REPLACEMENTS:
                parts.append(FIELD_LABEL_REPLACEMENTS[field])
            else:
                parts.append(field)

        string = '%s\n' % '\t'.join(parts)
        self.buffer.append(string)


    def next(self):
        with self.buffer_lock:
            try:
                line = self.buffer.pop(0)
            except IndexError:
                # empty buffer
                line = None

            # fill buffer
            if self.nEOF and len(self.buffer) < self.buffer_size and not self.buffer_thread:
                self.buffer_thread = BufferThread(self)
                self.buffer_thread.start()

        # wait for buffer if its completely empty
        if self.nEOF and not line:
            self.buffer_thread.join()
            with self.buffer_lock:
                line = self.buffer.pop(0)

        if line:
            return line

        raise StopIteration

    def __iter__(self):
        return self