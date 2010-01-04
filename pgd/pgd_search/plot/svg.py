""" Set of classes for helping deal with SVG graphics """


class SVG():
    """
    Class is a container for basic drawing objects.  It is used to record
    actions that make up an svg so that they may be redrawn using Jquery, pycairo
    or another drawing tool
    """

    def __init__(self):
        self.operations = []
        #self.rects = []
        #self.lines = []
        #self.texts = []

    def line(self,x,y,x1,y1,stroke=1,color='black'):
        self.operations.append(Line(x,y,x1,y1,stroke,color))

    def rect(self, x,y,height,width,stroke=1,color='black',fill=None,data=None):
        self.operations.append(Rect(x,y,height,width,stroke,color, fill, data))

    def text(self, x,y,text, size=16,fontfamily='Verdana', fill='black'):
        self.operations.append( Text(x,y,text,size,fontfamily, fill))

    def to_dict(self):
        """
        Convert object to dictionary so that it may be json serialized
        """
        return [op.__dict__ for op in self.operations]


class Line():
    def __init__(self, x,y,x1,y1,stroke=1,color='black'):
        self.type = 'line'
        self.x = x
        self.y = y
        self.x1 = x1
        self.y1 = y1
        self.stroke = stroke
        self.color = color


class Rect():
    def __init__(self,x,y,height,width,stroke=0,color='black',fill=None,data=None):
        self.type = 'rect'
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.stroke = stroke
        self.color = color
        self.fill = fill
        self.data = data


class Text():
    def __init__(self, x,y,text, size,color='#000000', fontfamily='Verdana'):
        self.type = 'text'
        self.x = x
        self.y = y
        self.text = text
        self.size = size
        self.fontfamily = fontfamily
        self.color = color