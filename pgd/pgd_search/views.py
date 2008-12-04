from django.http import HttpResponse, HttpResponseRedirect
from django.template import RequestContext, Context, loader
from django.conf import settings
from django.shortcuts import render_to_response
from django import forms
import math

from ConfDistFuncs import *
from svg import *

from constants import AA_CHOICES, SS_CHOICES
from pgd_search.models import *
from pgd_core.models import *


residueCountChoices = []
for i in range(1,searchSettings.segmentSize+1):
    residueCountChoices.append((i,i))

class SearchFormBase(forms.Form):
    threshold       = forms.ChoiceField(choices=[(25,25),(90,90)], required=False)
    resolutionMin   = forms.FloatField(required=False, min_value=0, initial=0, widget=forms.TextInput(attrs={'size':3}))
    resolutionMax   = forms.FloatField(required=False, min_value=0, initial=1.5, widget=forms.TextInput(attrs={'size':3}))
    proteins        = forms.ModelMultipleChoiceField(queryset=Protein.objects.all().order_by('code'), required=False)
    residues        = forms.ChoiceField(choices=residueCountChoices, initial=5)

# Build a dict for the fields of variable number
form_dict = {'__module__' : 'pgd_search.views'}

start = 0 - (searchSettings.segmentSize-1) / 2
stop  = int(math.ceil((searchSettings.segmentSize-1) / 2.0))+1
residueIndexes = range(start, stop, 1)
for i in residueIndexes:
    form_dict["aa_%i" % i]      = forms.ChoiceField(choices=AA_CHOICES, required=False, widget=forms.SelectMultiple(attrs={'class':'field'}))
    form_dict["aa_i_%i" % i]    = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

    form_dict["ss_%i" % i]      = forms.ChoiceField(choices=SS_CHOICES, required=False, widget=forms.Select(attrs={'class':'field'}))
    form_dict["ss_i_%i" % i]    = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

    # the loops here are just to save on space/typing
    for j in range(1,8):
        form_dict["a%i_%i" % (j,i)]     = forms.ChoiceField(choices=AA_CHOICES, required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["a%i_i_%i" % (j,i)]   = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))
    for j in range(1,6):
        form_dict["L%i_%i" % (j, i)]    = forms.FloatField(required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["L%i_i_%i" % (j, i)]  = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))
    for j in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta"):
        form_dict["%s_%i" % (j, i)]     = forms.FloatField(required=False, widget=forms.TextInput(attrs={'class':'field', 'size':8}))
        form_dict["%s_i_%i" % (j, i)]   = forms.IntegerField(required=False, widget=forms.HiddenInput(attrs={'class':'include'}))

# Create the Search Form with the fields from the dict
SearchForm = type('SearchForm', (SearchFormBase,), form_dict)

def search(request):
    if request.method == 'POST': # If the form has been submitted
        form = SearchForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            search = processSearchForm(form)
            return HttpResponseRedirect('/results/') # Redirect after POST
    else:
        form = SearchForm() # An unbound form

    #construct a list of values for i
    iValues = []
    for i in residueIndexes:
        if i < 0:
            iValues.append((i,'%i'%i))
        elif i == 0:
            iValues.append((i,'i'))
        else:
            iValues.append((i,'+%i'%i))

    #order the residue properties in way that django template can handle it better
    residueFields = []
    for i in residueIndexes:
        dict = {}
        for prefix in ("ss", "aa", "phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            dict[prefix] =  form['%s_%i' % (prefix, i)]
            dict['%s_i' % prefix] =  form['%s_i_%i' % (prefix, i)]
            dict['index'] = i
        residueFields.append(dict)

    return render_to_response('search.html', {
        'MEDIA_URL': settings.MEDIA_URL,
        'form': form,
        'maxLength' : searchSettings.segmentSize,
        'iValues':iValues,
        'residueFields':residueFields
    })


"""
Encode a list of AA choices into an integer value
"""
def encodeAA(list):
    aa = 0
    for i in range(len(AA_CHOICES)):
        if AA_CHOICES[i][0] in list:
            # bitwise or to add value
            aa = aa | (1 << i)
        i += 1
    return aa


"""
Decode an integer into a list of AA choices
"""
def decodeAA(val):
    list = []
    for i in range(23):
        # bitwise shift check value
        if (val & (1 << i)) != 0:
            list.append(AA_CHOICES[i][0])

    return list


"""
Process a search form copying its data into a search object
"""
def processSearchForm(form):
    data = form.cleaned_data
    length = searchSettings.segmentSize

    #create a new search object
    search = Search()

    #get protein properties
    search.residueCount  = int(data['residues'])
    search.resolutionMin = float(data['resolutionMin'])
    search.resolutionMax = float(data['resolutionMax'])

    #save search object so its residue parameters can be added
    search.save()

    #get list of proteins to filter
    for value in data['proteins']:
        searchCode = Search_code()
        searchCode.code = value
        search.codes.add(searchCode)

    #process per residue properties
    start = 0 - (search.residueCount-1) / 2
    stop  = int(math.ceil((search.residueCount-1) / 2.0))+1
    for i in range(start, stop, 1):
        residue = Search_residue()
        residue.index   = i

        #process ss
        residue.ss      = data['ss_%i' % i]

        #process aa
        residue.aa_int  = encodeAA([data['aa_%i' % i]])

        #process all other fields
        for prefix in ("phi", "psi", "ome", "chi", "bm", "bs", "bg", "h_bond_energy", "zeta", 'a1','a2','a3','a4','a5','a6','a7','L1','L2','L3','L4','L5'):
            residue.__dict__[prefix] = data['%s_%i' % (prefix, i)]
            residue.__dict__['%s_include' % prefix] = data['%s_i_%i' % (prefix, i)]

        search.residues.add(residue)


""" 
Renders a conformational distribution graph
@return: retusns an SVG instance.
"""
def drawGraph(xStart=-180, yStart=-180, xEnd=180, yEnd=180, attribute='Observations', xProperty='phi', yProperty='psi', reference=None, residue=0, xBin=10, yBin=10):
    svg = SVG()

    x = 55;
    y = 45;
    height = 400;
    width = 400;
    hashsize = 10

    #background
    svg.rect(x, y, width, height, 1, '#00fffff', '#222222');
    #svg.rect(0, 0, width+90, height+90, 1, '#00fffff');
    #svg.rect(x, y, width, height, 1, '#666666');

    #border
    svg.rect(x, y, width, height, 1, '#000000');

    #axis
    svg.line( x, y+height/2, x+width, y+height/2, 1, '#666666');
    svg.line( x+width/2, y, x+width/2, y+height, 1, '#666666');

    #hashes
    for i in range(9):
        hashx = x+(width/8.0)*i
        hashy = y+(height/8.0)*i
        svg.line( hashx, y+height, hashx, y+height+hashsize, 1, '#000000');
        svg.line( x, hashy, x-hashsize, hashy, 1, '#000000');

    #x axis text
    xtext  = -180
    xtext1 = 180
    step = (xtext1 - xtext) / 4
    for i in range(5):
        text = xtext + step*i
        hashx = x+(width/4)*i-(4*len(str(text)))
        svg.text(hashx, y+height+hashsize*3, str(text),12)

    #y axis text
    ytext  = -180
    ytext1 = 180
    step = (ytext1 - ytext) / 4
    for i in range(5):
        text = xtext + step*i
        hashy = y+(height/4)*i+7
        svg.text(x-20-(8*len(str(text))), hashy, str(text),12)

    #title text
    len1 = 220 - len(xProperty)*7/2 - len(xProperty)*7/2
    len2 = 182 - len(attribute)*7/2
    svg.text(len1,15, 'Plot of %s vs. %s' % (xProperty,yProperty), 12)
    svg.text(len2,35, 'Shading Based Off of %s' % attribute, 12)

#ob = 140
#


    cdp = ConfDistPlot(
            400,            #height
            400,            #width
            0,              #Xpadding
            0,              #Ypadding
            x,              #Xoffset
            y,              #Yoffset
            xStart,         #Xstart
            xEnd,           #Xend
            yStart,         #Ystart
            yEnd,           #Yend
            xBin,           #Xbin
            yBin,           #Ybin
            xProperty,      #X property
            yProperty,      #Y property
            '1sny',         #protein code
            attribute       #property
    )

    boxes = cdp.Plot()
    return (svg,boxes)

def RGBTuple(rgbString):
    sub = rgbString[-6:]
    red = int(sub[:2],16)/255.0
    green = int(sub[2:4],16)/255.0
    blue = int(sub[4:], 16)/255.0
    return (red,green,blue)

def line(input, context):
    context.move_to(input.x, input.y)
    context.line_to(input.x1, input.y1)
    context.set_line_width(input.stroke)
    r,g,b = RGBTuple(input.color)
    context.set_source_rgba(r,g,b,1)
    context.stroke()

def rect(input, context):
    context.rectangle(input.x, input.y, input.width, input.height)

    if input.fill:
        r,g,b = RGBTuple(input.fill)
        context.set_source_rgba(r,g,b,1)
        context.fill()

    if input.color:
        red, green, blue = RGBTuple(input.color)
        context.set_source_rgba(red,green,blue,1)
        context.set_line_width(input.stroke)
        context.stroke()

"""
render the conf dist graph to a png and return it as the response
this results in the image being downloaded by the user
"""
def renderToPNG(request):
    import cairo

    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, boxes = drawGraph(
                        data['x'],
                        data['y'],
                        data['x1'],
                        data['y1'],
                        data['attribute'],
                        data['xProperty'],
                        data['yProperty'],
                        data['reference'],
                        data['residue'],
                        data['xBin'],
                        data['yBin'])

    else:
        form = PlotForm() # An unbound form
        svg,boxes = drawGraph()

    width = 500
    height = 500

    response = HttpResponse(mimetype="image/png")
    response['Content-Disposition'] = 'attachment; filename="plot.png"'
    surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
    ctx = cairo.Context (surface)

    for rec in svg.rects:
        rect(rec, ctx)

    for action in svg.lines:
        line(action, ctx)

    surface.write_to_png(response)

    return response



ATTRIBUTE_CHOICES = [
                    ("Observations","Observations"),
                    ("L1","L1"),
                    ("L2","L2"),
                    ("L3","L3"),
                    ("L4","L4"),
                    ("L5","L5"),
                    ("a1","a1"),
                    ("a2","a2"),
                    ("a3","a3"),
                    ("a4","a4"),
                    ("a5","a5"),
                    ("a6","a6"),
                    ("a7","a7"),
                    ("ome","ome"),
                    ("chi","chi"),
                    ]
                    
PROPERTY_CHOICES = [
                    ("L1","L1"),
                    ("L2","L2"),
                    ("L3","L3"),
                    ("L4","L4"),
                    ("L5","L5"),
                    ("a1","a1"),
                    ("a2","a2"),
                    ("a3","a3"),
                    ("a4","a4"),
                    ("a5","a5"),
                    ("a6","a6"),
                    ("a7","a7"),
                    ("ome","ome"),
                    ("chi","chi"),
                    ("phi","phi"),
                    ("psi","psi"),
                    ]

class PlotForm(forms.Form):
    attribute       = forms.ChoiceField(choices=ATTRIBUTE_CHOICES, initial='Observations')
    xProperty       = forms.ChoiceField(choices=PROPERTY_CHOICES, initial='phi')
    yProperty       = forms.ChoiceField(choices=PROPERTY_CHOICES, initial='psi')
    reference       = forms.FloatField(required=False, widget=forms.TextInput(attrs={'size':8}))
    x               = forms.IntegerField(initial=-180, widget=forms.TextInput(attrs={'size':4}))
    x1              = forms.IntegerField(initial=180, widget=forms.TextInput(attrs={'size':4}))
    y               = forms.IntegerField(initial=-180, widget=forms.TextInput(attrs={'size':4}))
    y1              = forms.IntegerField(initial=180, widget=forms.TextInput(attrs={'size':4}))
    residue         = forms.ChoiceField(choices=[(i,'i') if i == 0 else (i,i) for i in range(start,stop)], initial=0)
    xBin            = forms.IntegerField(initial=10, widget=forms.TextInput(attrs={'size':4}))
    yBin            = forms.IntegerField(initial=10, widget=forms.TextInput(attrs={'size':4}))


"""
render conf dist plot using jquery.svg
"""
def renderToSVG(request):
    if request.method == 'POST': # If the form has been submitted
        form = PlotForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            data = form.cleaned_data
            svg, boxes = drawGraph(
                        data['x'],
                        data['y'],
                        data['x1'],
                        data['y1'],
                        data['attribute'],
                        data['xProperty'],
                        data['yProperty'],
                        data['reference'],
                        data['residue'],
                        data['xBin'],
                        data['yBin'])

    else:
        form = PlotForm() # An unbound form
        svg,boxes = drawGraph()

    return render_to_response('graph.html', {
        'MEDIA_URL': settings.MEDIA_URL,
        'form': form,
        'svg': svg,
        'boxes': boxes,
        'referenceValues' : RefDefaults()
    })
