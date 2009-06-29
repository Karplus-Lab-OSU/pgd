
   #--------------------------------------------------------------------------------------------------------------------
   # File: CleanDB.php
   # Purpose: Classes and defs associated with plotting data. 
   # Author: Mike Marr
   # Date: 9/28/05
   # Use: Use ConfDistPlot to create a plot, other classes are used by ConfDistPlot
   #-------------------------------------------------------------------------------------------------------------------

from pgd_core.models import *
from pgd_search.models import *
from constants import *
import math

ANGLES = ('ome', 'phi', 'psi', 'chi', 'zeta')

NON_FIELDS = ('Observations', 'all')

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
    def __init__(self, xSize, ySize, xPadding, yPadding, xOffset, yOffset, xMin, xMax, yMin, yMax, xbin, ybin, xText, yText, ref, residue, querySet):

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
        self.residue = int(math.ceil(searchSettings.segmentSize/2.0)-1) + (residue if residue else 0)

        # Printf-style string for the given residue
        self.resString = "r%i_%%s"%self.residue

        # Graphed quantity and its field string representation, e.g. 'r4_ome'
        self.ref = ref
        self.refString = self.resString%ref

        # Dictionary of bins, keyed by a tuple of x-y coordinates in field units
        #   i.e. (<x value>, <y value>)
        self.bins = {}

        # Labels for graph axes
        self.xText,self.yText = xText,yText
        self.xTextString,self.yTextString = self.resString%xText,self.resString%yText

        # Variable to store number of values in the bin with the most values
        self.maxObs = 0

        # Pick fields for retrieving values
        self.fields     = [(field,self.resString%str(field)) for field in (
            [
                xText, yText,
            ] if ref == "Observations" else [
                field for field,none in PLOT_PROPERTY_CHOICES
            ] if ref == "all" else [
                xText, yText, ref
            ]
        )]

        # 
        self.querySet = querySet.exclude(reduce(
            # exclude dummy values
            lambda x,y: x|y,
            [Q(**{"%s__in"%fieldString:(999.90,0)}) for field,fieldString in self.fields]
        )).filter(
            # Exclude values outside the plotted values
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
        xMarkOff,yMarkOff = int(xMin/xbin),int(yMin/ybin)

        for entry in self.querySet.values(*[fieldString for field,fieldString in self.fields]):
            
            # Adjustments for the axes values
            xAdj,yAdj = int(math.floor(entry[self.xTextString] / xbin)),int(math.floor(entry[self.yTextString] / ybin))

            key = (xAdj,yAdj)
            xDex,yDex = xAdj - xMarkOff,yAdj - yMarkOff
            if xText in ANGLES: xDex = xDex%xModder
            if yText in ANGLES: yDex = yDex%yModder

            if self.bins.has_key(key):
                # Append the observation to the bin entry...
                self.bins[key]['obs'].append(entry)
            else:
                # ...or add a new entry to the bins dict

                self.bins[key] = {
                    'obs'         : [entry],
                    'pixCoords'   : {
                        # The pixel coordinates of the x and y values
                        'x' : xMarks[xDex],
                        'y' : yMarks[yDex],
                        'width'  : widths[xDex],
                        'height' : heights[yDex],
                    }
                }

        # Calculate stats for each bin
        for bin in self.bins.values():

            obs = bin['obs']
            bin['count'] = len(obs)
            
            # Find the bin with the most observations
            if self.maxObs < bin['count']:
                self.maxObs = bin['count']

            # Calculate bin stats for each field
            for field,fieldString in self.fields:

                # Skip axes fields, unless 'll' is the indicated field
                #  (This reduces unnecessary stats calculations)
                if self.ref != 'all' and field in (self.xTextString, self.yTextString): continue

                # If only 1 observation, avg is the only value and stdev is zero
                if bin['count'] == 1:
                    bin['%s_avg'%fieldString] = obs[0][fieldString]
                    bin['%s_std'%fieldString] = 0

                # If >1 observation, calculate the avg and stdev
                else:
                    bin['%s_avg'%fieldString], bin['%s_std'%fieldString] = (
                        # Use the appropriate stats method
                        getCircularStats if field in ANGLES else getLinearStats
                    )(
                        [ob[fieldString] for ob in obs],
                        bin['count'],
                    )

    # ******************************************************
    # Plots observations
    # ******************************************************
    def Plot(self):
        
        binVals = []
        
        # Calculate stats regarding the distribution of averages in cells
        if self.ref not in NON_FIELDS and len(self.bins):
            if self.ref in ANGLES:
                meanPropAvg,stdPropAvg = getCircularStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                stdPropAvgX3 = 180 if stdPropAvg > 60 else 3*stdPropAvg
            else:
                meanPropAvg,stdPropAvg = getLinearStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                minPropAvg = meanPropAvg - 3*stdPropAvg
                maxPropAvg = meanPropAvg + 3*stdPropAvg

        # Color the bins
        for key in self.bins:

            bin = self.bins[key]
            num = bin['count']

            if self.ref in NON_FIELDS:
                scale = math.log(num+1, self.maxObs+1)
                color = map(
                    lambda x: x*scale,
                    (255.0,180.0,200.0)
                )
            elif self.ref in ANGLES:
                avg = bin['%s_avg'%self.refString]
                straight = avg - meanPropAvg
                difference = (
                    straight
                ) if -180 < straight < 180 else (
                    (360 if straight < 0 else -360) + straight
                )
                if -difference >= stdPropAvgX3 or difference >= stdPropAvgX3:
                    color = [255,-75,255]
                else:
                    scale = 0.5+((
                            math.log(
                                difference+1,
                                stdPropAvgX3+1
                            )
                        ) if difference >= 0 else (
                            -math.log(
                                -difference+1,
                                stdPropAvgX3+1
                          )
                       ))/2
                    color = map(
                        lambda x: x*scale,
                        (255.0,180.0,200.0)
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
                        (255.0,180.0,200.0)
                    )
            color[1] += 75

            #convert decimal RGB into HEX rgb
            fill = ''.join('%02x'%round(x) for x in color)

            # add rectangle to list
            binVals.append(
                [
                    bin['pixCoords']['x'],
                    bin['pixCoords']['y'] - bin['pixCoords']['height'],
                    bin['pixCoords']['width'],
                    bin['pixCoords']['height'],
                    fill,
                    fill,
                    bin,
                    key,
                ] + ([0,0] if self.ref in NON_FIELDS else [bin['%s_avg'%self.refString],bin['%s_std'%self.refString]])
            )

        return binVals


    # *****************************************************************************
    # Prints out the query results in a dump file
    # *****************************************************************************
    def PrintDump(self, writer):

        residue = self.residue

        #fields to include, order in this list is important
        STATS_FIELDS = ('phi','psi','ome','L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7','chi','zeta','h_bond_energy')
        avgString = '%s_avg'%self.resString
        stdString = '%s_avg'%self.resString
        STATS_FIELDS_STRINGS = reduce(
            lambda x,y:x+y,
            ((avgString%stat,stdString%stat) for stat in STATS_FIELDS)
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

        #output the generic titles
        for title in STATS_FIELDS:
            if title in field_replacements:
                title = field_replacements[title]
            elif not title in lower_case_fields:
                title = title.capitalize()
            writer.write('\t')
            writer.write('%sAvg' % title)
            writer.write('\t')
            writer.write('%sDev' % title)

        # Cycle through the binPoints
        for key in self.bins:
            bin = self.bins[key]
            writer.write('\n')

            # x axis range
            writer.write(key[0])
            writer.write('\t')
            writer.write(key[0]+self.xbin)

            # y-axis range
            writer.write('\t')
            writer.write(key[1])
            writer.write('\t')
            writer.write(key[1]+self.ybin)

            # observations
            writer.write('\t')
            writer.write(bin['count'])

            # Start averages and standard deviations
            for fieldStat in STATS_FIELDS_STRINGS:
                writer.write('\t')
                writer.write(round(bin[fieldStat], 1))


# ******************************************************
# Returns default reference values
# ******************************************************
def RefDefaults():
    return {
                'phi': {'ref':180, 'stepsize':1, 'custom':False},
                'Observations': [-10, -7, -4, -1, 2, 5, 8, 10, 10 ],
                'L1': { 
                        'ref':1.330, 
                        'stepsize':0.0025, 
                        'custom': 
                        False},
                'L2': { 
                        'ref': 1.465, 
                        'custom':False,
                        'stepsize':0.005},
                'L3': { 
                        'ref': 1.530,
                        'custom':False,
                        'stepsize':0.005},
                'L4': { 
                        'ref':1.525,
                        'custom':False,
                        'stepsize':0.005},
                'L5': { 
                        'ref': 1.240,
                        'custom':False,
                        'stepsize':0.005},
                'L6': { 
                        'ref': 1.330,
                        'custom':False,
                        'stepsize':0.005},
                'L7': { 
                        'ref': 1.465,
                        'custom':False,
                        'stepsize':0.005},
                'a1': { 
                        'ref': 121,
                        'custom':False,
                        'stepsize':1},
                'a2': { 
                        'ref': 110,
                        'custom':False,
                        'stepsize':1},
                'a3': { 
                        'ref': 110,
                        'custom':False,
                        'stepsize':1},
                'a4': { 
                        'ref': 110,
                        'custom':False,
                        'stepsize':1},
                'a5': { 
                        'ref': 120,
                        'custom':False,
                        'stepsize':0.5},
                'a6': { 
                        'ref': 117,
                        'custom':False,
                        'stepsize':1},
                'a7': { 
                        'ref': 123,
                        'custom':False,
                        'stepsize':1},
                'ome': { 
                        'ref': 180,
                        'custom':False,
                        'stepsize':1}
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
