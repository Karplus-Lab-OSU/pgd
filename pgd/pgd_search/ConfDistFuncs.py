
   #--------------------------------------------------------------------------------------------------------------------
   # File: CleanDB.php
   # Purpose: Classes and defs associated with plotting data. 
   # Author: Mike Marr
   # Date: 9/28/05
   # Use: Use ConfDistPlot to create a plot, other classes are used by ConfDistPlot
   #-------------------------------------------------------------------------------------------------------------------

from pgd_core.models import *
from constants import *
import math

#-------------------------------------------------------------------------------------------------------------------
# Stores a coordinate including actualy x,y values and x,y translated to pixel space
# Coordinates also have a dat string which is a unique identifier for the coordinate
class Coord():

    x   =  None # actual x,y values and pixel space x,y
    y   = None
    xP  = None
    yP  = None 
    dat = None                # ID and CODE forming a unique key for the object

    # ******************************************************
    # Constructor
    # ******************************************************
    def __init__(self, x, y, xP, yP, dat):
        self.dat = dat
        self.x = x
        self.y = y
        self.xP = xP
        self.yP = yP

#-------------------------------------------------------------------------------------------------------------------
# Class that divides an array of coordinates
# into bins specificed by an area
#-------------------------------------------------------------------------------------------------------------------
class Bin():

    # ******************************************************
    # Constructor
    # ******************************************************
    def __init__(self, xLen, yLen, xMin, xMax, yMin, yMax, points):
        self.xLen = xLen    # size of bounding box
        self.yLen = yLen
        self.xMin = xMin # Boundaries for x and y
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.bins   = {}        # Bins
        self.points = points # list of points
        self.maxObs = 0 # initial max observed
        self.MakeBins() # Create the bins

    # ******************************************************
    # Returns the bin coordinates for a specified location
    # ******************************************************
    def GetBinCoords(self, x, y):
        # Find the bin by dividing the coordinate by the bounding box length, rounding down
        # and then re-applying the length. I'm not sure why I did self ...
        xFloor = math.floor( x / self.xLen )
        xFloor = self.xLen * xFloor

        yFloor = floor ( y / self.yLen )
        yFloor = self.yLen * yFloor

        # Create a coordinate for the x and y bins that specificies the minimum and maximum
        # for the bin
        coords[0] = Coord(xFloor, xFloor + self.xLen, 0, 0, NULL) 
        coords[1] = Coord(yFloor, yFloor + self.yLen, 0, 0, NULL)

        return coords


    # ******************************************************
    # Returns the bin at a specificied coordinate
    # x and y represent x and y pixel location
    # ******************************************************
    def GetBin(self, x, y):
        xBin = math.floor( x / self.xLen )
        yBin = math.floor( y / self.yLen )
        key = '%s-%s' % (xBin,yBin)
        return self.bins[key]

    # ******************************************************
    # Get the number of observations at a particular bin
    # x and y represent a bin location
    # ******************************************************
    def GetObs(self, x, y):
        bin = self.GetBin(x, y)
        return bin.numObs

    # ******************************************************
    # Makes the bins
    # ******************************************************
    def MakeBins(self):
        # Iterate through all the points in the plot that need to be binned

        for i in range(len(self.points)):
            # Determine which bin a particular point belongs to
            xBin = math.floor( self.points[i].x / self.xLen )
            yBin = math.floor( self.points[i].y / self.yLen )

            key = '%s-%s' % (xBin,yBin)

            # If no bin currently exists, make a new binpoint

            if self.bins.has_key(key):
                # add to the bin a new observation
                self.bins[key].AddObs( self.points[i] )
            #otherwise create the bin
            else:
                self.bins[key] = BinPoint( self.points[i], self.points[i].xP, self.points[i].yP, xBin, yBin )

            # if the bin now holds the max number of observations, increased the max num of observations for ANY bin
            if self.bins[key].numObs > self.maxObs:
                self.maxObs = self.bins[key].numObs

#-------------------------------------------------------------------------------------------------------------------
# Class that contains information regarding an individual bin in a plot
#-------------------------------------------------------------------------------------------------------------------
class BinPoint():

    numObs    = None    # Number of observations
    avg       = None    # average for the bin    
    sum       = None    # sum of all data
    deviation = None    # standard deviation of data            
    obs       = []    # Array of observations
    xP        = None    # X in pixel space
    yP        = None    # Y in pixel space
    colorStep = None    # The color that should be used for the bin
    xBin      = None    # the bin's location in the grand scheme of things
    yBin      = None

    # ******************************************************
    # Saves original data x, y and pixel space x1 y1
    # and which xbin ybin the point belongs to
    # ******************************************************
    def __init__(self, coord, xP, yP, xBin, yBin):
        self.xP = xP
        self.yP = yP
        self.numObs = 0
        self.avg = None
        self.colorStep = None
        self.AddObs(coord)    # add the original x,y data to the list of observations for self bin
        self.xBin = xBin
        self.yBin = yBin


    # ******************************************************
    # Returns the x bin #
    # ******************************************************
    def GetXBin(self):
        return self.xBin


    # ******************************************************
    # Returns the y bin #
    # ******************************************************
    def GetYBin(self):
        return self.yBin


    # ******************************************************
    # Sets color bin
    # ******************************************************
    def SetColorStep(self, color):
        self.colorStep = color


    # ******************************************************
    # Gets color bin
    # ******************************************************
    def GetColor(self):
        return self.colorStep

    # ******************************************************
    # Saves original data x, y into the observation array
    # ******************************************************
    def AddObs(self, coord):
        self.obs.append(coord)
        self.numObs = self.numObs + 1


    # ******************************************************
    # Computes stats(Avg, standard deviation) for the bin and a specific stat, such as phi, psi, chi ...
    # ******************************************************
    def ComputeStats(self, key):
        # Average
        sum = 0
        for i in range(self.numObs):
            # By adding 360 to omega when it's less than -90, we shift half of the tall peak
            # over to the far right of a -180 to +180 omega plot, into the 180-270 range.
            if key == "ome":
                if self.obs[i].dat[key] < -90:
                    self.obs[i].dat[key] = ( self.obs[i].dat[key] + 360 )

            sum += self.obs[i].dat[key]

        self.avg[key] = sum / self.numObs

        # Std Deviation
        sum = 0
        for i in range(self.numObs):
            sum += pow( (self.obs[i].dat[key] - self.avg[key]), 2 )

        if self.numObs > 1:
            self.dev[key] = sqrt(sum / ( self.numObs - 1 ))

        

    # ******************************************************
    # Returns the average for a specificed value, such as phi, psi, L1 ...
    # ******************************************************
    def GetAvg(self, key):
        if self.avg[key] == None:
            self.ComputeStats(key)

        return self.avg[key]


    # ******************************************************
    # Returns the standard deviation for a specified value    
    # ******************************************************
    def GetDev(self, key):
        if self.dev[key] == None:
            self.ComputeStats(key)

        return self.dev[key]



#-------------------------------------------------------------------------------------------------------------------
# Class that plots conformation distribution plots
# Construction: X, Y image dimensions
#                             X, Y offsets from top left corner
#                             Query to used to populate plot
#-------------------------------------------------------------------------------------------------------------------
class ConfDistPlot():

    xSize       = None      # Size of image horizontally
    ySize       = None      # Size of image vertically
    xOff        = None      # Horizontal offset to start plot
    yOff        = None      # Vertical offset to start plot
    img         = None      # the image itself
    black       = None      # the color black
    white       = None      # the color white
    queryResult = None      # the sql query result
    theQuery    = None      # Query used
    points      = None      # all data points
    plotBin     = None      # bins for the plot
    ref         = None      # Reference attribute for shading

    xRange      = None      # X min and Max
    yRange      = None      # Y min and Max
    xPlotSize   = None      # X dimension for plot
    yPlotSize   = None      # Y dimension for plot
    xPixelSize  = None      # size of one x pixel
    yPixelSize  = None      # size of one y pixel

    xText       = None      # text for x axis
    yText       = None      # test for y axis
    fontSize    = None      # size of font
    fontHeight  = None      # height of font
    fontWidth   = None      # width of font

    maxObs      = None      # maximum number of observations
    REF         = None      # Reference values
    USEREF      = None
    table       = None      # table to use
    xbin        = None
    ybin        = None

    # ******************************************************
    # Constructor
    # x, y:             size of plot
    # x0, y0:         offset from top right corner of image for plot to being
    # Min, Max:    min and max values of plot
    # text:            x and y axis text
    # query:            sql query to use to get data for the plot
    # ref:                reference attribute for shading
    # ******************************************************
    def __init__(self, x, y, xPadding, yPadding, xOffset, yOffset, xMin, xMax, yMin, yMax, xbin, ybin, xText, yText, code, ref):
        self.xSize = x
        self.ySize = y
        self.xOff = xOffset
        self.yOff = yOffset
        self.xPadding = xPadding
        self.yPadding = yPadding
        self.ref = ref
        self.xbin = xbin
        self.ybin = ybin

        # Plotting region is the space allowed inside the image
        # excluding the border area of offsets
        self.xPlotSize = self.xSize - 2 * self.xPadding
        self.yPlotSize = self.ySize - 2 * self.yPadding
        self.xRange = [xMin, xMax]
        self.yRange = [yMin, yMax]

        self.xPixelSize = (self.xRange[1] - self.xRange[0]) / float(self.xPlotSize)    # Range / PlotSize = PixelSize
        self.yPixelSize = (self.yRange[1] - self.yRange[0]) / float(self.yPlotSize)
        self.xText = xText
        self.yText = yText
        self.maxObs = 0
        self.points = []
        self.bin = None
        self.plotBin = None

        #self.theQuery = query
        #self.SetQuery(query)        # Runs the query
        #self.MakeBasicImage()        #Creates the basic image
        #self.table = table
        self.code = code

        # Reference array that contains information for a specific type of value
        self.REF = {
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
        self.USEREF = self.REF

    # ******************************************************
    # Returns reference values used
    # key: value of interest
    # ******************************************************
    def GetColorRefs(self, key):
        arr[0] = self.USEREF[key]['ref']
        arr[1] = self.USEREF[key]['stepsize']
        return arr

    # ******************************************************
    # Sets reference value for a certain key
    # key: variable
    # attrib: attribute
    # value: value to set to
    # ******************************************************
    def SetRefAttrib(self, key, attrib, value):
        self.USEREF[key][attrib] = value
        self.USEREF[key]['custom'] = true


    # ******************************************************
    # Resets reference value for a certain key
    # key: variable
    # ******************************************************
    def ResetRefAttrib(self, key):
        self.USEREF[key] = self.REF[key]

    # ******************************************************
    # Creates a color array 
    # colorInterval: Interval for colors 
    # ******************************************************
    def CreateColors(self, colorInterval):

        COLORS = {
                        -10 :[colorInterval * 0,colorInterval * 0,colorInterval * 0],
                         -9 :[colorInterval * 1,colorInterval * 1,colorInterval * 1],
                         -8 :[colorInterval * 2,colorInterval * 2,colorInterval * 2],
                         -7 :[colorInterval * 3,colorInterval * 3,colorInterval * 3],
                         -6 :[colorInterval * 4,colorInterval * 4,colorInterval * 4],
                         -5 :[colorInterval * 5,colorInterval * 5,colorInterval * 5],
                         -4 :[colorInterval * 6,colorInterval * 6,colorInterval * 6],
                         -3 :[colorInterval * 7,colorInterval * 7,colorInterval * 7],
                         -2 :[colorInterval * 8,colorInterval * 8,colorInterval * 8],
                         -1 :[colorInterval * 9,colorInterval * 9,colorInterval * 9],
                         0  :[ colorInterval * 10,colorInterval * 10,colorInterval * 10],
                         1  :[ colorInterval * 11,colorInterval * 11,colorInterval * 11],
                         2  :[ colorInterval * 12,colorInterval * 12,colorInterval * 12],
                         3  :[ colorInterval * 13,colorInterval * 13,colorInterval * 13],
                         4  :[ colorInterval * 14,colorInterval * 14,colorInterval * 14],
                         5  :[ colorInterval * 15,colorInterval * 15,colorInterval * 15],
                         6  :[ colorInterval * 16,colorInterval * 16,colorInterval * 16],
                         7  :[ colorInterval * 17,colorInterval * 17,colorInterval * 17],
                         8  :[ colorInterval * 18,colorInterval * 18,colorInterval * 18],
                         9  :[ colorInterval * 19,colorInterval * 19,colorInterval * 19],
                         10 :[ colorInterval * 20,colorInterval * 20,colorInterval * 20]
                       }
        return COLORS


    # ******************************************************
    # Determines color coding for a bin
    # key: value of interest
    # val: value of the value
    # ******************************************************
    def DetermineColor(self, key, val):

        colorMax = 255 # Maximum color value
        colorInterval = round(255 / 21, 0 )    # intervals of color
        #COLORS = self.CreateColors(colorInterval)

        # use ranges for RGB to introduce colorful steps
        COLORS = {}
        redStart = 0;
        greenStart = 75;
        blueStart = 0;
        redInterval = round((255-redStart)/21)
        greenInterval = round((255-greenStart)/21)
        blueInterval = round((200-blueStart)/21)
        for i in range(21):
            COLORS[i-10] = [
                i*redInterval+redStart,
                i*greenInterval+greenStart,
                i*blueInterval+blueStart,
            ]

        colorStep = 0

        #Observations setup
        if key == 'Observations':
            if val == 1:
                    colorStep = self.REF['Observations'][1]
            elif val == 2:
                    colorStep = self.REF['Observations'][2]
            elif val == 5:
                    colorStep = self.REF['Observations'][3]
            elif val == 10:
                    colorStep = self.REF['Observations'][4]
            elif val == 25:
                    colorStep = self.REF['Observations'][5]
            elif val == 50:
                    colorStep = self.REF['Observations'][6]
            elif val == 100:
                    colorStep = self.REF['Observations'][7]
            elif val == 250:
                    colorStep = self.REF['Observations'][8]
            elif val == 500:
                    colorStep = self.REF['Observations'][9]
            elif val == 1000:
                    colorStep = self.REF['Observations'][10]
            else:
                    colorStep = 0

        else:
            variance = val - self.USEREF[key]['ref'] # Variance from the reference value
            colorStep = round(variance / self.USEREF[key]['stepsize'], 0 ) # color to be used based on how far from teh reference value
            if  colorStep < -10:
                colorStep = -10    # Bounds on the color
            if colorStep > 10:
                colorStep = 10

        color = COLORS[colorStep]    # Color to be used
        color.append(colorStep)        # Store color bin after the actual colors, I do self to save the color into the bin at a later point

        return color


    # ******************************************************
    # Plots observations
    # ******************************************************
    def PlotPoints(self, svg=None):
        if not svg:
            svg = SVG()

        # Only bins need to be painted, so cycle through the bins
        for key in self.plotBin.bins:
                bin = self.plotBin.bins[key]

                # Determine actual image location of the bin
                
                x = bin.xBin * self.plotBin.xLen
                y = bin.yBin * self.plotBin.yLen
                xC = ((x  - (self.xRange[0])) / self.xPixelSize + self.xOff)
                yC = ((-1 * y + self.yRange[0]) / self.yPixelSize + self.yPlotSize + self.yOff)

                
                if x < self.xRange[0] or x > self.xRange[1]:
                    continue
                if y < self.yRange[0] or y > self.yRange[1]:
                    continue

                # Bins are an area of pixel space, find the rectangle that describes
                # the area the bin uses
                xMin = math.floor(xC)
                xMax = xMin + round(self.plotBin.xLen / self.xPixelSize)
                yMax = math.floor(yC)
                yMin = yMax - round(self.plotBin.yLen / self.yPixelSize)

                # The number of observations in the bin and the max number of observations
                num = self.plotBin.GetObs(x, y)

                # As long as there is an observation, plot the rectangle for the bin
                if( num >= 1 ):
                    # Special Case for # obs
                    if self.ref == "Observations":
                        color = self.DetermineColor( self.ref, num)
                    else:
                        color = self.DetermineColor( self.ref, self.plotBin.bins[key].GetAvg(self.ref) )

                    self.plotBin.bins['%s-%s'%(bin.xBin,bin.yBin)].SetColorStep(color[-1])

                    #convert decimal RGB into HEX rgb
                    fill = '%s%s%s' % (hex(int(color[0]))[2:], hex(int(color[1]))[2:], hex(int(color[2]))[2:])

                    # add rectangle to svg
                    svg.rect( xMin+1, yMin+1, xMax-2-xMin, yMax-2-yMin, 1, color='#%s' % fill, fill='#%s' % fill )

        return svg


    # ******************************************************
    # Plots the points
    # ******************************************************
    def Plot(self, svg=None):

        # Turn all the query results into an array of points
        protein = Protein.objects.filter(code=self.code)[0]
        for residue in protein.residues.all():
        #while( row = mysql_fetch_array(self.queryResult) )

            #code = residue.code
            #id = residue.id

            xOrig = residue.__dict__[self.xText]    # Original Values of X and Y
            yOrig = residue.__dict__[self.yText]

            x = (xOrig  - (self.xRange[0])) / self.xPixelSize + self.xOff
            y = (yOrig + self.yRange[0]) / self.yPixelSize + self.yPlotSize + self.yOff

            self.points.append(Coord(xOrig, yOrig, x, y, residue.id))

        # Create a bins for the values
        self.plotBin = Bin(self.xbin, self.ybin, self.xRange[0], self.xRange[1], self.yRange[0], self.yRange[1], self.points)
        self.maxObs = self.plotBin.maxObs
        # Plot the bad boy
        return self.PlotPoints(svg)


    # *******************************************************************************
    # Prints out the query results in a dump file
    # *******************************************************************************
    '''def PrintDump(self):

        # Determine the file name for a dump
        handle = NULL
        out = NULL
        fileName = session_id() . "-plot.txt"

        echo "<center><p><a href=\"dump/fileName\">Download info</a></p></center>"
        handle = fopen("/var/www/pgd/database/dump/" . fileName, "w")

        out .= "PhiStart" . "\t" . "PhiStop" . "\t"
        out .= "PsiStart" . "\t" . "PsiStop" . "\t"
        out .= "Observations" . "\t"
        out .= "PhiAvg" . "\t"
        out .= "PhiDev" . "\t"
        out .= "PsiAvg" . "\t"
        out .= "PsiDev" . "\t"
        out .= "L1Avg" . "\t"
        out .= "L1Dev" . "\t"
        out .= "L2Avg" . "\t"
        out .= "L2Dev" . "\t"
        out .= "L3Avg" . "\t"
        out .= "L3Dev" . "\t"
        out .= "L4Avg" . "\t"
        out .= "L4Dev" . "\t"
        out .= "L5Avg" . "\t"
        out .= "L5Dev" . "\t"
        out .= "a1Avg" . "\t"
        out .= "a1Dev" . "\t"
        out .= "a2Avg" . "\t"
        out .= "a2Dev" . "\t"
        out .= "a3Avg" . "\t"
        out .= "a3Dev" . "\t"
        out .= "a4Avg" . "\t"
        out .= "a4Dev" . "\t"
        out .= "a5Avg" . "\t"
        out .= "a5Dev" . "\t"
        out .= "a6Avg" . "\t"
        out .= "a6Dev" . "\t"
        out .= "a7Avg" . "\t"
        out .= "a7Dev" . "\t"
        out .= "OmeAvg" . "\t"
        out .= "OmeDev" . "\t"
        out .= "\n"

        # Cycle through the bin areas
        # Derived from MakeMap()
        while( list(xBin, yArr) = each(self.plotBin.bins) )
            while( list(yBin, binVal) = each(yArr) )
                xLen = self.plotBin.xLen
                yLen = self.plotBin.yLen

                # x and y points
                x = xBin * xLen
                y = yBin * yLen

                if( x < self.xRange[0] || x > self.xRange[1] )print svg
                    break
                if( y < self.yRange[0] || y > self.yRange[1] )
                    break

                box = self.plotBin.GetBinCoords(x, y)

                out .= box[0].x . "\t" . box[0].y . "\t" # Phi range
                out .= box[1].x . "\t" . box[1].y . "\t" # Psi range
                out .= binVal.numObs . "\t" # Observations
                # Start averages and standard deviations
                out .= round(binVal.GetAvg('phi'), 1) . "\t" . round(binVal.GetDev('phi'), 3) . "\t" # Phi
                out .= round(binVal.GetAvg('psi'), 1) . "\t" . round(binVal.GetDev('psi'), 3) . "\t" # Psi
                out .= round(binVal.GetAvg('L1'), 3) . "\t" . round(binVal.GetDev('L1'), 3) . "\t" # L1
                out .= round(binVal.GetAvg('L2'), 3) . "\t" . round(binVal.GetDev('L2'), 3) . "\t" # L2
                out .= round(binVal.GetAvg('L3'), 3) . "\t" . round(binVal.GetDev('L3'), 3) . "\t" # L3
                out .= round(binVal.GetAvg('L4'), 3) . "\t" . round(binVal.GetDev('L4'), 3) . "\t" # L4
                out .= round(binVal.GetAvg('L5'), 3) . "\t" . round(binVal.GetDev('L5'), 3) . "\t" # L5
                out .= round(binVal.GetAvg('a1'), 3) . "\t" . round(binVal.GetDev('a1'), 3) . "\t" # a1
                out .= round(binVal.GetAvg('a2'), 3) . "\t" . round(binVal.GetDev('a2'), 3) . "\t" # a2
                out .= round(binVal.GetAvg('a3'), 3) . "\t" . round(binVal.GetDev('a3'), 3) . "\t" # a3
                out .= round(binVal.GetAvg('a4'), 3) . "\t" . round(binVal.GetDev('a4'), 3) . "\t" # a4
                out .= round(binVal.GetAvg('a5'), 3) . "\t" . round(binVal.GetDev('a5'), 3) . "\t" # a5
                out .= round(binVal.GetAvg('a6'), 3) . "\t" . round(binVal.GetDev('a6'), 3) . "\t" # a6
                out .= round(binVal.GetAvg('a7'), 3) . "\t" . round(binVal.GetDev('a7'), 3) . "\t" # a7
                out .= round(binVal.GetAvg('ome'), 1) . "\t" . round(binVal.GetDev('ome'), 3) . "\t" # Omega
                out .= "\n"

        fwrite( handle, out )
'''

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