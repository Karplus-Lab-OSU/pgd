
   #--------------------------------------------------------------------------------------------------------------------
   # File: CleanDB.php
   # Purpose: Classes and defs associated with plotting data. 
   # Author: Mike Marr
   # Date: 9/28/05
   # Use: Use ConfDistPlot to create a plot, other classes are used by ConfDistPlot
   #-------------------------------------------------------------------------------------------------------------------

import math

from django.db.models import Count, Avg, StdDev

from pgd_constants import *
from pgd_core.models import *
from pgd_search.models import *
from pgd_search.statistics.aggregates import DirectionalAvg, DirectionalStdDev

ANGLES = ('ome', 'phi', 'psi', 'chi', 'zeta')
NON_FIELDS = ('Observations', 'all')

"""
COLOR_RANGES - an RGB color setting for determining the range of colors in a plot
Made up of the MAX values for each Red, Green, Blue.  Plus an adjustment for each
RGB value.  A scale will be applied to each number equally then the adjustment added
This causes a range with all colors at the same proportion.  The adjustment causes
a grouping of colors closer to the max value.
"""
COLOR_RANGES = { 
    'green':(
        (255.0,180.0,200.0),
        (0,75,0)
     ),
    'blue':(
        (180.0,200.0,180.0),
        (0,0,75)
     ),
    'red':(
        (130.0, 200.0, 180.0),
        (115,0,0)
    ),
    'black':(
        (180.0,180.0,180.0),
        (75,75,75)
     )
}


# getCircularStats: returns (average, standard deviation) of a list of values
#   values: a list of the values to be examined
#   size:   the size of the list (in case it has been counted elsewhere)
def getCircularStats(values,size):

    # Store locals for speed
    lsin = math.sin
    lcos = math.cos
    lradians = math.radians
    lpow = math.pow

    # Circular Average - use some fancy trig that takes circular values
    #   into account.  This requires all values to be converted to radians.
    values = filter(lambda x:x!=None, values)
    size = len(values)

    if size == 1:
        return values[0],0

    radAngles = [lradians(val) for val in values]
    radAvg = math.atan2(
        sum([lsin(radAngle) for radAngle in radAngles])/size,
        sum([lcos(radAngle) for radAngle in radAngles])/size,
    )

    # Standard Deviation - shift the range of deviations +180 by applying
    #   %(2*pi) to all angles.  This creates a range of deviations -180-540.
    #   Values greater than 180 are then shifted back by substracting from
    #   360, resulting in deviations -180-180.  From there the Stdev formula
    #   is the same.
    msum = 0
    lpi = math.pi
    lpi_2 = lpi*2
    for radAngle in radAngles:
        straight = radAngle%lpi_2 - radAvg
        msum += lpow(straight if straight < lpi else lpi_2 - straight, 2)

    return math.degrees(radAvg),math.degrees(math.sqrt(msum/(size-1)))

# getLinearStats: returns (average, standard deviation) of a list of values
#   values: a list of the values to be examined
#   size:   the size of the list (in case it has been counted elsewhere)
def getLinearStats(values,size):

    # Average
    avg = sum(values)/size

    # Standard Deviation
    return avg,math.sqrt(
        sum([
            pow(value - avg, 2)
            for value in values
        ])/(size-1)
    )

#-------------------------------------------------------------------------------------------------------------------
# Class that plots conformation distribution plots
# Construction: X, Y image dimensions
#                             X, Y offsets from top left corner
#                             Query to used to populate plot
#-------------------------------------------------------------------------------------------------------------------
class ConfDistPlot():

    # ******************************************************
    # Constructor
    # Size:     size of plot
    # Padding:  space to either side of the plot
    # 0ffset:   offset from top right corner of image for plot to being
    # Min, Max: min and max field values of plot
    # bin:      field value size of a bin
    # Text:     axes of the plot
    # ref:      the field of interest (if 'all', observation count is plotted)
    # residue:  index of the residue of interest (-n...-1,0,1...n)
    # querySet: Django queryset
    # ******************************************************
    def __init__(self, xSize, ySize, xPadding, yPadding, xOffset, yOffset, xMin, xMax, yMin, yMax, xbin, ybin, xText, yText, ref, sigmaVal, residue_attribute, residue_xproperty, residue_yproperty, querySet, color='green'):
        self.color = color
        self.sigmaVal = sigmaVal

        # Convert unicode to strings
        xText,yText,ref = str(xText),str(yText),str(ref)

        # Width/height in field units of graph bins
        self.xbin = xbin
        self.ybin = ybin

        # Difference between the possible min and max axes values
        xLimit = xMax - xMin
        yLimit = yMax - yMin
        #   Adjustments for circular quantities
        if xText in ANGLES:
            xModder = int(360/xbin)
            if xMax < xMin:
                xLimit = xLimit%360
        if yText in ANGLES:
            yModder = int(360/ybin)
            if yMax < yMin:
                yLimit = yLimit%360

        # Width/height of a graph bin in pixels
        self.width  = round(self.xbin/((xLimit) / float(xSize - 2 * xPadding)))
        self.height = round(self.ybin/((yLimit) / float(ySize - 2 * yPadding)))

        # Index of the residue of interest in the segment
        self.residue_attribute = residue_attribute
        self.residue_xproperty = residue_xproperty
        self.residue_yproperty = residue_yproperty

        # get field prefix for this residue
        self.resString, self.refString = self.create_res_string(self.residue_attribute, ref)
        self.resXString, self.xTextString = self.create_res_string(self.residue_xproperty, xText)
        self.resYString, self.yTextString = self.create_res_string(self.residue_yproperty, yText)
        self.ref = ref

        # Dictionary of bins, keyed by a tuple of x-y coordinates in field units
        #   i.e. (<x value>, <y value>)
        self.bins = {}

        # Labels for graph axes
        self.xText,self.yText = xText,yText

        # Variable to store number of values in the bin with the most values
        self.maxObs = 0

        # index set creation
        self.index_set = set([self.resString,self.resXString,self.resYString])

        # Pick fields for retrieving values
        if ref == "Observations":
            self.fields = [(xText,self.xTextString), (yText,self.yTextString)]
            self.stats_fields = []
        elif ref == "all":
            self.fields = [(field,i%(str(field))) for field in ([field for field,none in PLOT_PROPERTY_CHOICES]) for i in self.index_set]
            self.stats_fields = self.fields
        else:
            self.fields = [(xText,self.xTextString), (yText,self.yTextString), (self.ref,self.refString)]
            self.stats_fields = [(self.ref,self.refString)]

        # Exclude values outside the plotted values
        self.querySet = querySet.filter(
            (Q(**{
                '%s__gte'%self.xTextString: xMin,
                '%s__lt'%self.xTextString: xMax,
            }) if (xMin <= xMax) else ( # Altered logic for circular values
                Q(**{'%s__gte'%self.xTextString: xMin}) |
                Q(**{'%s__lt'%self.xTextString: xMax})
            )) & (Q(**{
                '%s__gte'%self.yTextString: yMin,
                '%s__lt'%self.yTextString: yMax,
            }) if (yMin <= yMax) else ( # altered logic for circular values
                Q(**{'%s__gte'%self.yTextString: yMin}) |
                Q(**{'%s__lt'%self.yTextString: yMax})
            ))
        )

        # Total # of observations
        self.numObs = self.querySet.count()

        # create set of annotations to include in the query
        annotations = {'count':Count('id')}
        torsion_avgs = {}
        for field in self.fields:
            avg = '%s_avg' % field[1]
            stddev = '%s_stddev' % field[1]
            if field[0] in ANGLES:
                annotations[avg] = DirectionalAvg(field[1])
                torsion_avgs[field[0]] = {}
            else:
                annotations[avg] = Avg(field[1])
                annotations[stddev] = StdDev(field[1])
        annotated_query = self.querySet.annotate(**annotations)

        # determine aliases used for the table joins.  This is needed because
        # the aliases will be different depending on what fields were queried
        # even if the query is length 10, not all residues will be joined unless
        # each residue has a property in the where clause.
        x_alias = self.determine_alias(annotated_query, residue_xproperty)
        y_alias = self.determine_alias(annotated_query, residue_yproperty)
        attr_alias = self.determine_alias(annotated_query, residue_attribute)
        x_field = '%s.%s' % (x_alias, self.xText)
        y_field = '%s.%s' % (y_alias, self.yText)

        # calculating x,y bin numbers for every row.  This allows us
        # to group on the bin numbers automagically sorting them into bins
        # and applying the aggregate functions on them.
        x_aggregate = 'FLOOR(%s/%s)' % (x_field, xbin)
        y_aggregate = 'FLOOR(%s/%s)' % (y_field, ybin)
        annotated_query = annotated_query.extra(select={'x':x_aggregate, 'y':y_aggregate}).order_by('x','y')

        # add all the names of the aggregates and x,y properties to the list 
        # of fields to display.  This is required for the annotation to be
        # applied with a group_by.
        values = annotations.keys() + ['x','y']
        annotated_query = annotated_query.values(*values)

        # XXX remove the id field from the group_by.  By default django 
        # adds this to the group by clause.  This would prevent grouping
        # because the id is a unique field.  There is no official API for 
        # modifying group by and this is a big hack, but its a very simple
        # way of making this work
        annotated_query.query.group_by = []

        ### Sort the values from observations into bins
        # save these beforehand to avoid recalculating per bin
        xScale = xLimit / float(xSize - 2 * xPadding)
        yScale = yLimit / float(ySize - 2 * yPadding)
        yAllOffset = (ySize - 2 * yPadding) + yOffset
        #  Calculate the bin boundaries (to avoid recalculation)
        xMarks = [math.floor((mark*xbin) / xScale + xOffset)+1 for mark in range(0,int(math.floor(xLimit/xbin))+1)]
        widths = [xMarks[i+1] - xMarks[i] - 2 for i in range(0,len(xMarks)-1)]
        yMarks = [math.floor(-(mark*ybin) / yScale + yAllOffset)-1 for mark in range(0,int(math.floor(yLimit/ybin))+1)]
        heights = [yMarks[i] - yMarks[i+1] - 2 for i in range(0,len(yMarks)-1)]
        #  Adjust to make + indices for the Marks lists
        xMarkOff,yMarkOff = int(math.floor(xMin/xbin)), int(math.floor(yMin/ybin))
    
        for entry in annotated_query:

            x = int(entry['x'])
            y = int(entry['y'])
            key = (x,y)
            xDex = x - xMarkOff
            yDex = y - yMarkOff

            if xText in ANGLES: xDex = xDex%xModder
            if yText in ANGLES: yDex = yDex%yModder

            # add  entry to the bins dict
            bin = {
                'count' : entry['count'],
                'obs'         : [entry],
                'pixCoords'   : {
                    # The pixel coordinates of the x and y values
                    'x' : xMarks[xDex],
                    'y' : yMarks[yDex],
                    'width'  : widths[xDex],
                    'height' : heights[yDex],
                }
            }

            # add all statistics
            for k, v in entry.items():
                if k in ('x','y','count'):
                    continue
                bin[k] = v 

            if bin['count'] > 1:
                # if this is an angle the stddev must be calculated in separate query
                # using the circular standard deviation method
                #
                # due to null causing errors, each field must be run separate to filter
                # out nulls for just that field
                for field in self.stats_fields:
                    if field[0] in ANGLES:
                        if avg:
                            torsion_avgs[field[0]]["'%s:%s'" % (x,y)] = bin[avg]

            else:
                # no need for calculation, stddev infered from bincount
                for field in self.fields:
                    if field[0] in ANGLES:
                        bin['%s_stddev' % field[1]] = 0

            self.bins[key] = bin

            # Find the bin with the most observations
            if self.maxObs < bin['count']:
                self.maxObs = bin['count']

        # run queries to get stddevs for all torsion angles
        # this is done outside the main loop because all of the averages 
        # must first be collected so that they can be run in a single query
        #
        # This query uses a large case statement to select the average matching
        # the bin the result is grouped into.  This isn't the cleanest way
        # of doing this but it does work
        for field in self.stats_fields:
            if field[0] in ANGLES:
                stddev = '%s_stddev' % field[1]
                cases = ' '.join(['WHEN %s THEN %s' % (k,v) if v else '' for k,v in torsion_avgs[field[0]].items()])
                avgs = "CASE CONCAT(FLOOR(%s/10.0),':',FLOOR(%s/10.0)) %s END" % (x_field, y_field, cases)
                annotations = {stddev:DirectionalStdDev(field[1], avg=avgs)}
                bin_where_clause = ['NOT %s.%s IS NULL' % (attr_alias, field[0])]
                stddev_query = self.querySet \
                                    .extra(select={'x':x_aggregate, 'y':y_aggregate}) \
                                    .extra(where=bin_where_clause) \
                                    .annotate(**annotations) \
                                    .values(*annotations.keys()+['x','y']) \
                                    .order_by('x','y')
                stddev_query.query.group_by = []

                for r in stddev_query:
                    value = r[stddev] if r[stddev] else 0
                    self.bins[(r['x'],r['y'])][stddev] = value


    def create_res_string(self, index, property):
        """
        helper function for creating property references
        """
        if index == 0:
            prefix = ''
        elif index < 0:
            prefix = ''.join(['prev__' for i in range(index, 0)])
        else:
            prefix = ''.join(['next__' for i in range(index)])
        resString = '%s%%s' % prefix
        refString = '%s%s' % (prefix, property)
        
        return resString, refString
    
    
    def determine_alias(self, query, index):
        """
        determines the table alias used for a given residue index.
        
        XXX This takes into account django internal structure as of 12/29/2009
        this may change with future releases.
        
        query.join_map is a dict mapping a tuple of (table1, table2, fk, key)
        mapped to a list of aliases the table is joined on.  multiple aliases
        means the table was joined on itself multiple times.
        
        we must walk the list of joins to find the index number we want.
        
        @returns alias if table is joined, otherwise None
        """
        query = query.query
        if index == 0:
            return 'pgd_core_residue'
        if index > 0:
            k = ('pgd_core_residue','pgd_core_residue','next_id','id')
        else:
            k = ('pgd_core_residue','pgd_core_residue','prev_id','id')
            
        if not query.join_map.has_key(k):
            return None
        try:
            return query.join_map[k][int(math.fabs(index))-1]
        except IndexError:
            return None

        
            
    # ******************************************************
    # Plots observations
    # ******************************************************
    def Plot(self):
        binVals = []
        sig = self.sigmaVal
        # Calculate stats regarding the distribution of averages in cells
        if self.ref not in NON_FIELDS and len(self.bins):
            if self.ref in ANGLES:
                meanPropAvg,stdPropAvg = getCircularStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                stdPropAvgXSigma = 180 if stdPropAvg > 60 else sig*stdPropAvg
            else:
                meanPropAvg,stdPropAvg = getLinearStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                minPropAvg = meanPropAvg - sig*stdPropAvg
                maxPropAvg = meanPropAvg + sig*stdPropAvg

        colors, adjust = COLOR_RANGES[self.color]
        # Color the bins
        for key in self.bins:
            bin = self.bins[key]
            num = bin['count']

            if self.ref in NON_FIELDS:
                scale = math.log(num+1, self.maxObs+1)
                color = map(
                    lambda x: x*scale,
                    colors
                )
            elif self.ref in ANGLES:
                avg = bin['%s_avg'%self.refString]
                if avg and '%s_stddev'%self.refString in bin:
                    straight = avg - meanPropAvg
                    difference = (
                        straight
                    ) if -180 < straight < 180 else (
                        (360 if straight < 0 else -360) + straight
                    )
                else:
                    # no average, mark bin as outlier
                    difference = 9999

                if -difference >= stdPropAvgXSigma or difference >= stdPropAvgXSigma:
                    color = [255,-75,255]
                else:
                    scale = 0.5+((
                            math.log(
                                difference+1,
                                stdPropAvgXSigma+1
                            )
                        ) if difference >= 0 else (
                            -math.log(
                                -difference+1,
                                stdPropAvgXSigma+1
                          )
                       ))/2
                    color = map(
                        lambda x: x*scale,
                        colors
                    )
            else:
                avg = bin['%s_avg'%self.refString]
                if avg <= minPropAvg or avg >= maxPropAvg:
                    color = [255,-75,255]
                else:
                    scale = 0.5+((
                            math.log(
                                avg-meanPropAvg+1,
                                maxPropAvg-meanPropAvg+1
                            )
                        ) if avg > meanPropAvg else (
                            -math.log(
                                meanPropAvg-avg+1,
                                meanPropAvg-minPropAvg+1
                            )
                        ))/2
                    color = map(
                        lambda x: x*scale,
                        colors
                    )

            color[0] += adjust[0]
            color[1] += adjust[1]
            color[2] += adjust[2]

            #convert decimal RGB into HEX rgb
            fill = ''.join('%02x'%round(x) for x in color)

            # add rectangle to list
            if self.ref in NON_FIELDS:
                bin_avg, bin_stddev = 0,0 
            else:
                try:
                    bin_avg = bin['%s_avg'%self.refString]
                    bin_stddev = bin['%s_stddev'%self.refString]
                    if bin_avg == None:
                        continue
                except KeyError:
                    continue

            binVals.append(
                [
                    bin['pixCoords']['x'],
                    bin['pixCoords']['y'] - bin['pixCoords']['height'],
                    bin['pixCoords']['width'],
                    bin['pixCoords']['height'],
                    fill,
                    fill,
                    bin['count'],
                    key,
                    bin_avg,
                    bin_stddev
                ]
            )

        return binVals


    # *****************************************************************************
    # Prints out the query results in a dump file
    # *****************************************************************************
    def PrintDump(self, writer):

        residue = self.residue_attribute
        residueX = self.residue_xproperty
        residueY = self.residue_yproperty

        #fields to include, order in this list is important
        STATS_FIELDS = ('phi','psi','ome','L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7','chi','zeta')
        avgString = 'r%i_%s_avg'
        stdString = 'r%i_%s_stddev'

        index_set = set([residue,residueX,residueY])

        STATS_FIELDS_STRINGS = reduce(
            lambda x,y:x+y,
            ((avgString%(i,stat),stdString%(i,stat)) for stat in STATS_FIELDS for i in index_set),
        )

        lower_case_fields = ['a1','a2','a3','a4','a5','a6','a7']
        field_replacements = {
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
            'h_bond_energy':'HBond'
        }


        #capitalize parameters where needed
        if self.xText in lower_case_fields:
            xText = self.xText
        else:
            xText = self.xText.capitalize()

        if self.yText in lower_case_fields:
            yText = self.yText
        else:
            yText = self.yText.capitalize()

        #output the dynamic titles
        writer.write('%sStart' % xText)
        writer.write('\t')
        writer.write('%sStop' % xText)
        writer.write('\t')
        writer.write('%sStart' % yText)
        writer.write('\t')
        writer.write('%sStop' % yText)
        writer.write('\t')
        writer.write('Observations')

        residue_replacements = {
            0:u'(i-4)',
            1:u'(i-3)',
            2:u'(i-2)',
            3:u'(i-1)',
            4:u'(i)',
            5:u'(i+1)',
            6:u'(i+2)',
            7:u'(i+3)',
            8:u'(i+4)',
            9:u'(i+5)'
        }

        #output the generic titles
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            writer.write('\t')
            writer.write('%sAvg%s' % (title,residue_replacements[self.residue_attribute]))
            writer.write('\t')
            writer.write('%sDev%s' % (title,residue_replacements[self.residue_attribute]))

        #output the generic titles for x res
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            if len(index_set) > 1:
                writer.write('\t')
                writer.write('%sAvg%s' % (title,residue_replacements[self.residue_xproperty]))
                writer.write('\t')
                writer.write('%sDev%s' % (title,residue_replacements[self.residue_xproperty]))

        #output the generic titles for y res
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            if len(index_set) == 3:
                writer.write('\t')
                writer.write('%sAvg%s' % (title,residue_replacements[self.residue_yproperty]))
                writer.write('\t')
                writer.write('%sDev%s' % (title,residue_replacements[self.residue_yproperty]))

        # Cycle through the binPoints
        xbin = self.xbin
        ybin = self.ybin
        for key in self.bins:
            bin = self.bins[key]
            writer.write('\n')

            # x axis range
            writer.write(key[0]*xbin)
            writer.write('\t')
            writer.write((key[0]+1)*xbin)

            # y-axis range
            writer.write('\t')
            writer.write(key[1]*ybin)
            writer.write('\t')
            writer.write((key[1]+1)*ybin)

            # observations
            writer.write('\t')
            writer.write(bin['count'])

            # Start averages and standard deviations
            for fieldStat in STATS_FIELDS_STRINGS:
                writer.write('\t')
                val = bin[fieldStat] if fieldStat in bin else 0
                writer.write(val if val else 0)


# ******************************************************
# Returns default reference values
# ******************************************************
def RefDefaults():
    return {
                'phi': {
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'L1': {
                        'stepsize':'',
                        'min':'',
                        'max':'',},
                'L2': {
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'L3': {
                        'stepsize':'',
                        'min':'',
                        'max':''
                        },
                'L4': {
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'L5': {
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'L6': {
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'L7': {
                        'ref': 1.465,
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'a1': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a2': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a3': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a4': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a5': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a6': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'a7': {
                        'min':'',
                        'max':'',
                        'stepsize':''},
                'ome':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'chi':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'zeta':{
                        'min':-180,
                        'max':180,
                        'stepsize':10}
                }

if __name__ == "__main__":
    cdp = ConfDistPlot(400, 400, 100, 100,
                    -180, 180,
                    -180, 180,
                    10, 10,
                    "phi", "psi",
                    '1sny',
                    'Observations')

    svg = cdp.Plot()
    print svg.rects
