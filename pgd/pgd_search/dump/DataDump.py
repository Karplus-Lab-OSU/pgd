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
FIELDS = ['aa','a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5','ss','phi', 'psi', 'ome', 'chi1','chi2','chi3','chi4', 'bm', 'bs', 'bg', 'h_bond_energy', 'zeta']
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
    'a7':u'O-C-N(+1)',
    'a1_include':u'C(-1)-N-CA include',
    'a2_include':u'N-CA-CB include',
    'a3_include':u'N-CA-C include',
    'a4_include':u'CB-CA-C include',
    'a5_include':u'CA-C-O include',
    'a6_include':u'CA-C-N(+1) include',
    'a7_include':u'O-C-N(+1) include',
    'L1_include':u'C(-1)N include',
    'L2_include':u'N-CA include',
    'L3_include':u'CA-CB include',
    'L4_include':u'CA-C include',
    'L5_include':'C-O include',
    'phi_include':'phi include',
    'ome_include':'ome include',
    'chi1_include':'chi(1) include',
    'chi2_include':'chi(2) include',
    'chi3_include':'chi(3) include',
    'chi4_include':'chi(4) include',
    'bm_include':'bm include',
    'bs_include':'bs include',
    'bg_include':'bg include',
    'h_bond_energy_include':'H bond energy include',
    'zeta_include':'zeta include'
    }
FIELD_VALUE_REPLACEMENTS = {'aa':AA_CHOICES}
RESIDUE_FIELDS =    ['index','chainID','a1','a1_include','a2','a2_include','a3','a3_include','a4','a4_include','a5','a5_include','a6','a6_include','a7',
                    'a7_include','L1','L1_include','L2','L2_include','L3','L3_include','L4','L4_include','L5','L5_include','ss','phi', 'psi','phi_include', 'ome',
                    'ome_include','chi1','chi1_include','chi2','chi2_include','chi3','chi3_include','chi4','chi4_include', 'bm','bm_include','bs','bs_include','bg','bg_include', 'h_bond_energy','h_bond_energy_include', 'zeta','zeta_include']
SS_KEY_LIST = ['&alpha; helix','3<sub>10</sub> helix','&beta; sheet','Turn','Bend','&beta;-bridge','&pi; helix']
SS_HEADER = [u'Alpha Helix',u'3_10 Helix',u'Beta Sheet',u'Turn',u'Bend','Beta-Bridge','Pi Helix']

class BufferThread(Thread):
    """
    Adds processed segments to a Dump's buffer.  This class
    is used so that the buffering happens at the same time as the file
    being read.

    Page sizes are defined by the Dump class.  Each page contains a set
    of segments which are processed into text lines and added to the buffer
    """

    def __init__(self, parent):
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
                residue = segment
                if offset < 0:
                    while offset != 0:
                        residue = residue.prev
                        offset += 1
                elif offset > 0:
                    while offset != 0:
                        residue = residue.next
                        offset -= 1
                parts = [
                    str(self.parent.count) if first else ' ',
                    segment.protein_id,
                    string,
                    residue.oldID,
                    segment.chainID,
                ]
                #field values
                for field in FIELDS:
                    # replace field with display value if needed
                    if field in FIELD_VALUE_REPLACEMENTS:
                        code = residue.__dict__[field]
                        if code:
                            for k,v in FIELD_VALUE_REPLACEMENTS[field]:
                                if k == code:
                                    parts.append(str(v))
                    # just write value
                    else:
                        parts.append(str(residue.__dict__[field]))

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
        self.search = search
        self.pages = Paginator(search.querySet(), self.buffer_increment)
        self.page_max = self.pages.page_range[-1]
        self.current_page = 1
        self.count = 0
        self.nEOF = True
        self.create_meta_data(search)
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


    def create_meta_data(self,search):
        """Adds all of the relevent data about how the search was conducted."""
        #Add meta data begin tag to make parsing dumped searches easier
        self.buffer.append("***BEGIN_META_DATA***\n")
        #The first run sets up the headers
        parts = ['Dataset Date']
        for header in RESIDUE_FIELDS:
            if header in FIELD_LABEL_REPLACEMENTS:
                parts.append(str(FIELD_LABEL_REPLACEMENTS[header]))
            else:
                if header is 'ss':
                    parts+=SS_HEADER
                else:
                    parts.append(header)
        string = '%s\n' % '\t'.join(parts)
        self.buffer.append(string)
        parts = []
        #The rest of the loops fill in the data
        search_fields = search.residues.all()
        for residue in search_fields:
            parts.append(str(search.dataset_version))
            for key in RESIDUE_FIELDS:
                if key is 'ss':
                    ss = residue.__dict__[key]
                    for ss_key in SS_KEY_LIST:
                        parts.append(str(ss[ss_key]))
                else:
                    parts.append(str(residue.__dict__[key]))
            #At the end of a row, so join string and add newline
            string = '%s\n' % '\t'.join(parts)
            self.buffer.append(string)
            parts = []
        self.buffer.append("***END_META_DATA***\n")
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
        #print 'Here\'s the string!: '+string
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
