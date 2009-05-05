
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

# Takes angles between -180 and 180 and returns the
# angle of the shortest arc between the two
def shortCircle(first,second):
    straight = second - first
    return  (
                straight
            ) if -180 < straight < 180 else (
                (360 if straight < 0 else -360) + straight
            )

def getCircularStats(values,size):
    print "circular"

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


def getLinearStats(values,size):
    print "linear"

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
    # x, y:             size of plot
    # x0, y0:         offset from top right corner of image for plot to being
    # Min, Max:    min and max values of plot
    # text:            x and y axis text
    # query:            sql query to use to get data for the plot
    # ref:                reference attribute for shading
    # ******************************************************
    def __init__(self, xSize, ySize, xPadding, yPadding, xOffset, yOffset, xMin, xMax, yMin, yMax, xbin, ybin, xText, yText, ref, residue, querySet):

        #<john>
        xText,yText,ref = [str(field) for field in (xText,yText,ref)]
        self.xText = xText
        self.yText = yText
        self.ybin = ybin
        self.xbin = xbin
        self.yPixelSize = (yMax - yMin) / float(ySize - 2 * yPadding)
        self.xPixelSize = (xMax - xMin) / float(ySize - 2 * yPadding)
        self.residue = int(math.ceil(searchSettings.segmentSize/2.0)-1) + (residue if residue else 0)
        self.resString   = "r%i_%%s"%self.residue
        self.refString = self.resString%ref
        self.ref = ref
        self.bins   = {}
        self.xProperty,self.yProperty = self.resString%xText,self.resString%yText
        self.maxObs = 0

        # <firstloop>

        # pick fields to query
        self.fields     = [(field,self.resString%str(field)) for field in (
            [
                xText, yText,
            ] if ref == "Observations" else [
                field for field,none in PLOT_PROPERTY_CHOICES
            ] if ref == "all" else [
                xText, yText, ref
            ]
        )]
        self.querySet = querySet.exclude(reduce(
            lambda x,y: x|y,
            [Q(**{"%s__in"%fieldString:(999.90,0)}) for field,fieldString in self.fields]
        )).filter(
            Q(**{
                '%s__gte'%self.xProperty: xMin,
                '%s__lte'%self.xProperty: xMax,
                '%s__gte'%self.yProperty: yMin,
                '%s__lte'%self.yProperty: yMax,
            }) if (xMin <= xMax) else (
                (
                    (
                        (
                            Q(**{'%s__gte'%self.xProperty: xMin}) |
                            Q(**{'%s__lte'%self.xProperty: xMax})
                        ) & (
                            Q(**{'%s__gte'%self.yProperty: yMin}) |
                            Q(**{'%s__lte'%self.yProperty: yMax})
                        )
                    )
                )
            )
        )

        self.numObs = self.querySet.count()
        
        # save these beforehand to avoid recalculating per bin
        xScale = (xMax - xMin) / float(xSize - 2 * xPadding)
        yScale = (yMax - yMin) / float(ySize - 2 * yPadding)
        yAllOffset = (ySize - 2 * yPadding) + yOffset
        for entry in self.querySet.values(*[fieldString for field,fieldString in self.fields]):
            
            xAdj,yAdj = math.floor(entry[self.xProperty] / xbin),math.floor(entry[self.yProperty] / ybin)

            key = (xAdj,yAdj)

            if self.bins.has_key(key):
                self.bins[key]['obs'].append(entry)
            else:
                self.bins[key] = {
                    'obs'         : [entry],
                    'pixCoords'   : {
                        'x' :(xAdj*xbin - xMin) / xScale + xOffset,
                        'y' :(yMin - yAdj*ybin) / yScale + yAllOffset,
                    }
                }

        # </firstloop>
        # <secondloop>
        for bin in self.bins.values():
            obs = bin['obs']
            bin['count'] = len(obs)
            
            if self.maxObs < bin['count']:
                self.maxObs = bin['count']

            for field,fieldString in self.fields:
                if self.ref != 'all' and field in (self.xProperty, self.yProperty): continue
                # store functions and variables locally for speed optimization

                # if there is only 1 observation then avg is the only value and deviation is zero
                if bin['count'] == 1:
                    bin['%s_avg'%fieldString] = obs[0][fieldString]
                    bin['%s_std'%fieldString] = 0

                #if there is more than one observation then we must calculate the values
                else:
                    bin['%s_avg'%fieldString], bin['%s_std'%fieldString] = (
                        getCircularStats if field in (
                            'ome', 'phi', 'psi', 'chi', 'zeta'
                        ) else getLinearStats
                    )(
                        [ob[fieldString] for ob in obs],
                        bin['count'],
                    )
        # </secondloop>
        #</john>

    # ******************************************************
    # Plots observations
    # ******************************************************
    def Plot(self):
        
        binVals = []
        width   = round(self.xbin / self.xPixelSize)
        height  = round(self.ybin / self.yPixelSize)

        if width > 2:
            width -= 2
        else:
            width = 1

        if height > 2:
            height -= 2
        else:
            height = 1

        # Calculate stats regarding the distribution of averages in cells
        if self.ref not in ('Observations', 'all') and len(self.bins):
            if self.ref in ('ome', 'phi', 'psi', 'chi', 'zeta',):
                meanPropAvg,stdPropAvg = getCircularStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                stdPropAvgX3 = 180 if stdPropAvg > 60 else 3*stdPropAvg
            else:
                print self.refString,self.ref
                meanPropAvg,stdPropAvg = getLinearStats([bin['%s_avg'%self.refString] for bin in self.bins.values()], len(self.bins))
                minPropAvg = meanPropAvg - 3*stdPropAvg
                maxPropAvg = meanPropAvg + 3*stdPropAvg

        # Only bins need to be painted, so cycle through the bins
        for key in self.bins:
            bin = self.bins[key]

            # The number of observations in the bin and the max number of observations
            num = bin['count']

            # As long as there is an observation, plot the rectangle for the bin
            if num:

                if self.ref in ('Observations', 'all'):
                    scale = math.log(num+1, self.maxObs+1)
                    color = map(
                        lambda x: x*scale,
                        (255.0,180.0,200.0)
                    )
                elif self.ref in ('ome', 'phi', 'psi', 'chi', 'zeta',):
                    avg = bin['%s_avg'%self.refString]
                    difference = shortCircle(avg,meanPropAvg)
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
                # adjust dimensions to create 1 pixel border around bins
                # do not adjust if creating the border will result in a height less than 1

                # add rectangle to list
                binVals.append(
                    [
                        math.floor(bin['pixCoords']['x']) + (width > 2),
                        math.floor(bin['pixCoords']['y']) - round(self.ybin / self.yPixelSize) + (height > 2),
                        width,
                        height,
                        fill,
                        fill,
                        bin,
                        key,
                    ] + ([0,0] if self.ref in ('Observations', 'all') else [bin['%s_avg'%self.refString],bin['%s_std'%self.refString]])
                )

        return binVals


    # *******************************************************************************
    # Prints out the query results in a dump file
    # *******************************************************************************
    def PrintDump(self, writer):

        residue = self.residue

        #fields to include, order in this list is important
        STATS_FIELDS = ['phi','psi','ome','L1','L2','L3','L4','L5','a1','a2','a3','a4','a5','a6','a7','chi','zeta','h_bond_energy']
        lower_case_fields = ['a1','a2','a3','a4','a5','a6','a7']
        field_replacements = {'h_bond_energy':'HBond'}

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
            for field in STATS_FIELDS:
                writer.write('\t')
                writer.write(round(bin['%s_avg'%self.resString%field], 1))
                writer.write('\t')
                writer.write(round(bin['%s_std'%self.resString%field], 3))


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
