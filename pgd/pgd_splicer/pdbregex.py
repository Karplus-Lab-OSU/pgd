#!/usr/bin/env python

import re

lines = ['   90   1YKGA   146 -1.00  0.00     N   146   146      0     67     30  s34ul-fite (redu_ctase) [nadph] flavoprotein alpha- component'
         ,'   90   1YKGA   146 -1.25  0.00     N   146   146      0     67     30  sulfite (redu_ctase) [nadph] flavoprotein alpha- component'
         ,'   25   1CHC_    68 -1.00  0.00     N    68    68      0      8      8  Equ43ine he3rpes virus-1 (c3hc4, or ring domain) (NMR, 1 structure)'
         ]


#regex = '([-]?\d+\.\d+)'

#regex = '\s+(\d+)\s+([a-zA-Z0-9_])\s+(\d+)   \s+()   \s+()   \s+([a-zA-Z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+([a-zA-Z0-9_/(/) ]+)'

#WORKS
#regex = '\s+(\d+)\s+([a-zA-Z0-9_]+)\s+(\d+)\s+([-]?\d+\.\d+)\s+(\d+\.\d+)\s+([a-zA-Z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([a-zA-Z0-9/-_/(/) ]+)'

regex = '\s+(\d+)\s+([\w]+)\s+(\d+)\s+([-]?\d+\.\d+)\s+(\d+\.\d+)\s+([a-zA-Z]?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([,\w\-/(/) \[\]]+)'





p = re.compile(regex)
print p


for line in lines:
    m = p.match(line)
    print m.groups()