
   #--------------------------------------------------------------------------------------------------------------------
   # File: CleanDB.php
   # Purpose: Classes and defs associated with plotting data. 
   # Author: Mike Marr
   # Date: 9/28/05
   # Use: Use ConfDistPlot to create a plot, other classes are used by ConfDistPlot
   #-------------------------------------------------------------------------------------------------------------------

import math

import cairo
from django.db.models import Count, Avg, StdDev

from pgd_constants import *
from pgd_core.models import *
from pgd_search.models import *
from pgd_search.statistics.aggregates import DirectionalAvg, DirectionalStdDev, BinSort
from pgd_splicer.sidechain import sidechain_length_relationship_list, sidechain_angle_relationship_list
from svg import *

ANGLES = ('ome', 'phi', 'psi', 'chi1','chi2','chi3','chi4', 'zeta')
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


ROMAN_TO_GREEK = {
    'A':u'\u1D45',
    'B':u'\u1D5D',
    'G':u'\u03b3',
    'D':u'\u03b4',
    'E':u'\u03b5',
    #'H':'e',
    'Z':u'\u03B6',
    #'1':u'\u00B9',
    #'2':'\u00B2',
    #'3':'\u00B3',
    }


def format_atom(value):
    """ formats an atom bond (length or angle) using unicode characters """
    parts = []
    for i in value.split('_'):
        subscript = []
        for c in i[1:]:
            subscript.append(ROMAN_TO_GREEK[c] if c in ROMAN_TO_GREEK else c)
        parts.append('%s%s' % (i[0], ''.join(subscript)))
    return '-'.join(parts)


LABEL_REPLACEMENTS = {
            "L1":u'C\u207B\u00B9N',
            "L2":u'NC\u1D45',
            "L3":u'C\u1D45C\u1D5D',
            "L4":u'C\u1D45C',
            "L5":u'CO',
            "a1":u'C\u207B\u00B9NC\u1D45',
            "a2":u'NC\u1D45C\u1D5D',
            "a3":u'NC\u1D45C',
            "a4":u'C\u1D5DC\u1D45C',
            "a5":u'C\u1D45CO',
            "a6":u'C\u1D45CN\u207A\u00B9',
            "a7":u'OCN\u207A\u00B9',
            "ome":u'\u03C9',
            "chi1":u'\u03C7\u00B9',
            "chi2":u'\u03C7\u00B2',
            "chi3":u'\u03C7\u00B3',
            "chi4":u'\u03C7\u2084',
            "phi":u'\u03D5',
            "psi":u'\u03A8',
            'zeta':u'\u03B6',
            'h_bond_energy':'H Bond'
            }

for field in sidechain_angle_relationship_list:
    LABEL_REPLACEMENTS['sidechain_%s'%field] = '%s:%s' % (field[:3], format_atom(field[5:]))
    
for field in sidechain_length_relationship_list:
    LABEL_REPLACEMENTS['sidechain_%s'%field] = '%s:%s' % (field[:3], format_atom(field[5:]))


def getCircularStats(values,size):
    """
    getCircularStats: returns (average, standard deviation) of a list of values
    @param values: a list of the values to be examined
    @param size:   the size of the list (in case it has been counted elsewhere)
    """
    
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


def getLinearStats(values,size):
    """
    getLinearStats: returns (average, standard deviation) of a list of values
    @param values: a list of the values to be examined
    @param size:   the size of the list (in case it has been counted elsewhere)
    """

    # Average
    values = filter(lambda x: x!=None, values)
    avg = sum(values)/size

    # Standard Deviation
    return avg,math.sqrt(
        sum([
            pow(value - avg, 2)
            for value in values
        ])/(size-1)
    )


class ConfDistPlot():
    """ 
     Class that plots conformation distribution plots
     Construction: X, Y image dimensions
                                 X, Y offsets from top left corner
                                 Query to used to populate plot
    """ 

    
    def __init__(self, xSize, ySize, xMin, xMax, yMin, yMax,
                 xbin, ybin, xText, yText, ref, sigmaVal, residue_attribute,
                 residue_xproperty, residue_yproperty, querySet,
                 color='green',
                 background_color='#ffffff',
                 graph_color='#222222',
                 text_color='#000000',
                 hash_color='#666666'
                 ):
        """
         Constructor
         Size:     size of plot
         Padding:  space to either side of the plot
         0ffset:   offset from top right corner of image for plot to begining
         Min, Max: min and max field values of plot
         bin:      field value size of a bin
         Text:     field name to plot on each axis
         ref:      field name of interest (if 'all', observation count is plotted)
         residue:  index of the residue of interest (-n...-1,0,1...n)
         querySet: Django queryset
         
         color:    Hue to use for bin colorations
         background_color: color used for background of entire image
         graph_color: color used for background of plotted area
         text_color: color used for axis labels, hash labels, and title
         hash_color: color used for axis and hashes
        """
    
        # Convert unicode to strings
        xText,yText,ref = str(xText),str(yText),str(ref)
    
        # save properties
        self.querySet = querySet
        self.ref = ref
        self.xText = xText
        self.yText = yText
        self.x = xMin
        self.x1 = xMax
        self.y = yMin
        self.y1 = yMax
        self.sigmaVal = sigmaVal
        
        self.width = xSize
        self.height = ySize
        self.color = color
        self.background_color = background_color
        self.graph_color = graph_color
        self.text_color = text_color
        self.hash_color = hash_color    

        # Width/height in field units of graph bins
        self.xbin = xbin
        self.ybin = ybin

        # Difference between the possible min and max axes values
        xLimit = xMax - xMin
        yLimit = yMax - yMin
        # Index of the residue of interest in the segment
        self.residue_attribute = residue_attribute
        self.residue_xproperty = residue_xproperty
        self.residue_yproperty = residue_yproperty

    def query_bins(self):
        """
        Runs the query to calculate the bins and their relevent data
        """
        # local vars
        x = self.x
        x1 = self.x1
        y = self.y
        y1 = self.y1
        xbin = self.xbin
        ybin = self.ybin

        # if this is a circular spanning 180/-180 we must adjust the bins
        # so that the range falls between -180 and 180
        xlinear = True
        if x1 > 180:
            x = x1-360
            x1 = 180-self.x
            xlinear = False
        elif x1 < 0 and x > 0:
            x = x1
            x1 = self.x
            xlinear = False

        ylinear = True
        if y1 > 180:
            y = y1-360
            y1 = 180-self.y
            ylinear = False
        elif y1 < 0 and y > 0:
            y = y1
            y1 = self.y
            ylinear = False
            
        # Dictionary of bins, keyed by a tuple of x-y coordinates in field units
        #   i.e. (<x value>, <y value>)
        self.bins = {}

        # Variable to store number of values in the bin with the most values
        self.maxObs = 0

        # get field prefix for this residue
        self.resString, self.refString = self.create_res_string(self.residue_attribute, self.ref)
        self.resXString, self.xTextString = self.create_res_string(self.residue_xproperty, self.xText)
        self.resYString, self.yTextString = self.create_res_string(self.residue_yproperty, self.yText)

        # Exclude values outside the plotted values
        querySet = self.querySet.filter(
            (Q(**{
                '%s__gte'%self.xTextString: x,
                '%s__lte'%self.xTextString: x1,
            }) if (xlinear) else ( # Altered logic for circular values
                Q(**{'%s__gte'%self.xTextString: self.x}) |
                Q(**{'%s__lte'%self.xTextString: x})
            )) & (Q(**{
                '%s__gte'%self.yTextString: y,
                '%s__lte'%self.yTextString: y1,
            }) if (ylinear) else ( # altered logic for circular values
                Q(**{'%s__gte'%self.yTextString: self.y}) |
                Q(**{'%s__lte'%self.yTextString: y})
            ))
        )
        # Total # of observations
        self.numObs = querySet.count()
        
        # index set creation
        self.index_set = set([self.resString,self.resXString,self.resYString])
        
        # Pick fields for retrieving values
        if self.ref == "Observations":
            self.fields = [(self.xText,self.xTextString), (self.yText,self.yTextString)]
            self.stats_fields = []
        elif self.ref == "all":
            self.fields = [(field,i%(str(field))) for field in ([field for field,none in PLOT_PROPERTY_CHOICES]) for i in self.index_set]
            self.stats_fields = self.fields
        else:
            self.fields = [(self.xText,self.xTextString), (self.yText,self.yTextString), (self.ref,self.refString)]
            self.stats_fields = [(self.ref,self.refString)]

        # create set of annotations to include in the query
        annotations = {'count':Count('id')}
        torsion_avgs = {}
        for field in self.stats_fields:
            avg = '%s_avg' % field[1]
            stddev = '%s_stddev' % field[1]
            if field[0] in ANGLES:
                annotations[avg] = DirectionalAvg(field[1])
                torsion_avgs[field[0]] = {}
            else:
                annotations[avg] = Avg(field[1])
                annotations[stddev] = StdDev(field[1])
        annotated_query = querySet.annotate(**annotations)
        
        # sort and group by bins using an aggregate function that calculates
        # bin index based on bin size (in field units ie. degrees) and bin count.
        #
        # XXX django won't add aggregates to group by so to force it add the
        # create the aggregate functions and add them to a copy of the queryset
        # afterwards get the SQL from the aggregate function and add it as an
        # extra select clauses.  the select claus fields are added to the group by
        # statement.
        sortx = BinSort(self.xTextString, offset=x, bincount=xbin, max=x1)
        sorty = BinSort(self.yTextString, offset=y, bincount=ybin, max=y1)
        annotated_query.annotate(x=sortx, y=sorty)
        annotated_query = annotated_query.extra(select={'x':sortx.aggregate.as_sql(), 'y':sorty.aggregate.as_sql()})
        annotated_query = annotated_query.order_by('x','y')
        
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
        for entry in annotated_query:
            key = (int(entry['x']), int(entry['y']))
            # add  entry to the bins dict
            bin = {
                'count' : entry['count'],
                'obs'         : [entry],
                'pixCoords'   : key
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
                        if avg and bin[avg]:
                            torsion_avgs[field[0]]["'%s:%s'" % key] = bin[avg]

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
        prefix = self.get_prefix(self.residue_attribute)
        for field in self.stats_fields:
            
            if field[0] in ANGLES:
                stddev = '%s_stddev' % field[1]
                cases = ' '.join(['WHEN %s THEN %s' % (k,v) for k,v in filter(lambda x:x[1], torsion_avgs[field[0]].items())])
                if cases:
                    avgs = "CASE CONCAT(%s,':',%s) %s END" % (sortx.aggregate.as_sql(), sorty.aggregate.as_sql(), cases)
                    annotations = {stddev:DirectionalStdDev(field[1], avg=avgs)}
                    stddev_query = querySet \
                                    .extra(select={'x':sortx.aggregate.as_sql(), 'y':sorty.aggregate.as_sql()}) \
                                    .filter(**{'%s%s__isnull'%(prefix, field[1]):False}) \
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
        prefix = self.get_prefix(index)
        resString = '%s%%s' % prefix
        refString = '%s%s' % (prefix, property)
        return resString, refString

    def get_prefix(self, index):
        """ Helper for generating residue prefixes based on index """
        if index == 0:
            prefix = ''
        elif index < 0:
            prefix = ''.join(['prev__' for i in range(index, 0)])
        else:
            prefix = ''.join(['next__' for i in range(index)])
        return prefix

    def Plot(self):
        """
        Calculates and renders the plot
        """
        
        #cache local variables
        x = self.x
        y = self.y
        x1 = self.x1
        y1 = self.y1
        xbin = self.xbin
        ybin = self.ybin
        xText = self.xText
        yText = self.yText
        height = self.height
        width = self.width
        bg_color = self.background_color
        hash_color = self.hash_color
        text_color = self.text_color

        svg = SVG()

        # draw background
        #size ratio (470 = 1)
        ratio = width/560.0
    
        # offsets setup to give 415px for the graph for a default width of 520
        graph_x = round(width*.17857);
        graph_y = round(height*.11702);
        graph_height = height-2*graph_y;
        graph_width = width-2*graph_x;
        hashsize = 10*ratio
        
        # calculate bin count and sizes.
        if x>0 and x1<0:
            # account for circular coordinates
            xBinCount = math.ceil((360.0+x1-x)/xbin)
        else:
            xBinCount = math.ceil((float(x1)-x)/xbin)
        if y>0 and y1<0:
            # account for circular coordinates
            yBinCount = math.ceil((360.0+y1-y)/ybin)
        else:
            yBinCount = math.ceil((float(y1)-y)/ybin)
            
        # determine graph and bin sizes.  bins should always fall within the
        # borders of the graph, and the graph should always exactly fit the
        # bins.  to make this work the graph height/width must be adjusted to
        # fit the bins.  The unused portion will also be calculated
        binWidth = math.floor((graph_width-xBinCount+1)/xBinCount)
        binHeight = math.floor((graph_height-yBinCount+1)/yBinCount)
        graph_height_used = (binHeight+1)*yBinCount
        graph_width_used = (binWidth+1)*xBinCount
        unused = graph_height - graph_height_used;
        
        #image background
        svg.rect(0, 0, height+30, width, 0, bg_color, bg_color);
        #graph background
        svg.rect(graph_x, graph_y+unused, graph_height_used, graph_width_used, 0, self.graph_color, self.graph_color);
        #border
        svg.rect(graph_x+0.5, graph_y+0.5+unused, graph_height_used, graph_width_used, 1, hash_color);

        #draw data area (bins)
        self.query_bins()
        self.render_bins(svg, graph_x, graph_height+graph_y, binWidth, binHeight)

        #y axis
        if x < 0 and x1 > 0:
            xZero = (graph_width_used/(x1-x)) * abs (x)
            svg.line( graph_x+xZero, graph_y, graph_x+xZero, graph_y+graph_height_used, 1, hash_color);
        elif x > x1 :
            xZero = (graph_width_used/(360-abs(x1)-x)) * (180-x)
            svg.line( graph_x+xZero, graph_y, graph_x+xZero, graph_y+graph_height_used, 1, hash_color);
        #x axis
        if y < 0 and y1 > 0:
            yZero = graph_height_used+graph_y - (graph_height_used/(y1-y)) * abs (y)
            svg.line( graph_x, yZero+unused, graph_x+graph_width_used, yZero+unused, 1, hash_color);
        elif y > y1:
            yZero = graph_height_used+graph_y - (graph_height_used/(360-abs(y1)-y)) * (180-y)
            svg.line( graph_x, yZero+unused, graph_x+graph_width_used, yZero+unused, 1, hash_color);

        #hashes
        for i in range(9):
            hashx = graph_x+(graph_width_used/8.0)*i
            hashy = graph_y+(graph_height_used/8.0)*i
            svg.line( hashx, graph_y+graph_height, hashx, graph_y+graph_height+hashsize, 1, hash_color);
            svg.line( graph_x, hashy+unused, graph_x-hashsize, hashy+unused, 1, self.hash_color);
    
        #create a cairo surface to calculate text sizes
        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context (surface)
        ctx.set_font_size (12);
    
        #hash labels
        xstep = ((x1 - x)%360 if xText in ANGLES else (x1 - x))/ 4
        if not xstep: xstep = 90
        #ystep = (self.y1 - self.y) / 4
        ystep = ((y1 - y)%360 if yText in ANGLES else (y1 - y))/ 4
        if not ystep: ystep = 90
    
        #get Y coordinate for xaxis hashes, this is the same for all x-labels
        xlabel_y = graph_y+graph_height+hashsize*2.5
        ctx.set_font_size (12*ratio);
        for i in range(5):
            #text value
            xtext = ((x + xstep*i + 180)%360 - 180) if xText in ANGLES and x1 <= 180 else (x + xstep*i)
            #drop decimal if value is an integer
            xtext = '%i' % int(xtext) if not xtext%1 else '%.1f' %  xtext
            #get X coordinate of hash, offsetting for length of text
            xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(xtext)
            xlabel_x = graph_x+(graph_width_used/4)*i-xbearing-twidth
            #create label
            svg.text(xlabel_x, xlabel_y, xtext,12*ratio, text_color)
    
            #text value
            #ytext = self.y1 - ystep*i
            ytext = ((y + ystep*i + 180)%360 - 180) if yText in ANGLES and y1 <= 180 else (y + ystep*i)
            #drop decimal if value is an integer
            ytext = '%i' % int(ytext) if not ytext%1 else '%.1f' % ytext
            #get Y coordinate offsetting for height of text
            xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(ytext)
            ylabel_y = graph_y+(graph_height_used/4)*(4-i)-ybearing/2-theight/2+unused
            #Get X coordinate offsetting for length of hash and length of text
            ylabel_x = (graph_x-(ratio*28))-xbearing-twidth
            #create label
            svg.text(ylabel_x, ylabel_y, ytext,12*ratio, text_color)

        #title text
        xTitle = LABEL_REPLACEMENTS[xText] if xText in LABEL_REPLACEMENTS else xText
        yTitle = LABEL_REPLACEMENTS[yText] if yText in LABEL_REPLACEMENTS else yText
        title = 'Plot of %s vs. %s' % (xTitle,yTitle)
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(title)
        title_x = graph_x + (graph_width_used/2) - xbearing/2 - twidth/2
        svg.text(title_x,15*ratio, title, 12*ratio, text_color)
    
        attribute_title = LABEL_REPLACEMENTS[self.ref] if self.ref in LABEL_REPLACEMENTS else self.ref
        title = 'Shading Based Off of %s' % attribute_title
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(title)
        title_x = graph_x+(graph_width_used/2) - xbearing/2 - twidth/2
        svg.text(title_x,35*ratio, title, 12*ratio, text_color)
    
        #axis labels
        ctx.set_font_size (18*ratio);
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(xTitle)
        title_x = graph_x+(graph_width_used/2) - xbearing - twidth/2
        svg.text(title_x,graph_y+graph_height+hashsize*5, xTitle, 18*ratio, text_color)
    
        xbearing, ybearing, twidth, theight, xadvance, yadvance = ctx.text_extents(yTitle)
        title_y = (graph_x-(ratio*35))-xbearing-twidth+unused
        svg.text(25,graph_y+(graph_height/2)+xbearing+twidth/2, yTitle, 18*ratio, text_color, rotate=-90)

        return svg


    def render_bins(self, svg, xOffset, yOffset, binWidth, binHeight):
        """
        Renders the already calculated bins.
        """
        #cache variables
        sig = self.sigmaVal

        # Calculate stats regarding the distribution of averages in cells
        if self.ref not in NON_FIELDS and len(self.bins):
            key = '%s_avg'%self.refString
            values = [bin[key] for bin in filter(lambda x:x[key], self.bins.values())]
            if len(values):
                if self.ref in ANGLES:
                    meanPropAvg,stdPropAvg = getCircularStats(values, len(values))
                    stdPropAvgXSigma = 180 if stdPropAvg > 60 else sig*stdPropAvg
                else:
                    meanPropAvg,stdPropAvg = getLinearStats(values, len(values))
                    minPropAvg = meanPropAvg - sig*stdPropAvg
                    maxPropAvg = meanPropAvg + sig*stdPropAvg
            else:
                # no values, sigma is meaningless but set a value anyways the
                # remainder of the code will run.
                stdPropAvgXSigma = 0
                

        colors, adjust = COLOR_RANGES[self.color]
        # Color the bins
        for key in self.bins:
            bin = self.bins[key]
            num = bin['count']

            if self.ref in NON_FIELDS:
                scale = math.log(num+1, self.maxObs+1)
                color = map(
                    lambda z: z*scale,
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
            fill = '#%s' % ''.join('%02x'%round(x) for x in color)

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

            # create rect object.  positions are based on bin_index
            svg.rect(
                    bin['pixCoords'][0]*(binWidth+1)+xOffset+1,
                    yOffset-(bin['pixCoords'][1]+1)*(binHeight+1)+1,
                    binHeight,
                    binWidth,
                    0,
                    fill,
                    fill,
                    data = [
                        bin['count'],
                        key,
                        bin_avg,
                        bin_stddev
                    ]
            )


    def PrintDump(self, writer):
        """
        Prints out the query results in a dump file
        
        @param writer - any object that has a write(str) method
        """
        if not self.bins:
            self.query_bins()

        residue = self.residue_attribute
        residueX = self.residue_xproperty
        residueY = self.residue_yproperty

        #fields to include, order in this list is important
        STATS_FIELDS = ('phi','psi','ome','L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7','chi1','chi2','chi3','chi4','zeta')
        avgString = '%s%s_avg'
        stdString = '%s%s_stddev'
        print residue, residueX, residueY
        index_set = set([self.get_prefix(residue),self.get_prefix(residueX),self.get_prefix(residueY)])

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
        
        # convert index to label
        residue_attribute = '' if self.residue_attribute==0 else '%+d'%self.residue_attribute
        residue_xproperty = '' if self.residue_xproperty==0 else '%+d'%self.residue_xproperty
        residue_yproperty = '' if self.residue_yproperty==0 else '%+d'%self.residue_yproperty
        
        #output the generic titles
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            writer.write('\t')
            writer.write('%sAvg(i%s)' % (title,residue_attribute))
            writer.write('\t')
            writer.write('%sDev(i%s)' % (title,residue_attribute))

        #output the generic titles for x res
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            if len(index_set) > 1:
                writer.write('\t')
                writer.write('%sAvg(i%s)' % (title,residue_xproperty))
                writer.write('\t')
                writer.write('%sDev(i%s)' % (title,residue_xproperty))

        #output the generic titles for y res
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            if len(index_set) == 3:
                writer.write('\t')
                writer.write('%sAvg(i%s)' % (title,residue_yproperty))
                writer.write('\t')
                writer.write('%sDev(i%s)' % (title,residue_yproperty))

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


def RefDefaults():
    """
    Returns dictionary of default values for all properties.   These are used
    to provide defaults to fields that we do not want to automatically
    calculate dimensions for
    """
    return {
                'phi': {
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'L7': {
                        'ref': 1.465,
                        'stepsize':'',
                        'min':'',
                        'max':''},
                'ome':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'chi1':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'chi2':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'chi3':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'chi4':{
                        'min':-180,
                        'max':180,
                        'stepsize':10},
                'psi':{
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
