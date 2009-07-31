""" Set of classes for helping deal with SVG graphics """


"""
Class is a container for basic drawing objects.  It is used to record
actions that make up an svg so that they may be redrawn using Jquery, pycairo
or another drawing tool
"""
class SVG():
    def __init__(self):
        self.rects = []
        self.lines = []
        self.texts = []

    def line(self,x,y,x1,y1,stroke=1,color='black'):
        self.lines.append(Line(x,y,x1,y1,stroke,color))

    def rect(self, x,y,height, width, stroke=1,color='black', fill=None):
        self.rects.append(Rect(x,y,height, width, stroke,color, fill))

    def text(self, x,y,text, size=16,fontfamily='Verdana', fill='black'):
        self.texts.append( Text(x,y,text,size,fontfamily, fill))

    def to_dict(self):
        """
        Convert object to dictionary so that it may be json serialized
        """
        return {
            'lines': [line.__dict__ for line in self.lines],
            'rects': [rect.__dict__ for rect in self.rects],
            'texts': [text.__dict__ for text in self.texts]
        }

class Line():
    def __init__(self, x,y,x1,y1,stroke=1,color='black'):
        self.x = x
        self.y = y
        self.x1 = x1
        self.y1 = y1
        self.stroke = stroke
        self.color = color


class Rect():
    def __init__(self, x,y,height, width, stroke=1,color='black', fill=None):
        self.x = x
        self.y = y
        self.height = height
        self.width = width
        self.stroke = stroke
        self.color = color
        self.fill = fill


class Text():
    def __init__(self, x,y,text, size,color='#000000', fontfamily='Verdana'):
        self.x = x
        self.y = y
        self.text = text
        self.size = size
        self.fontfamily = fontfamily
        self.color = color
