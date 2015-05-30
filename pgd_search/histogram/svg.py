import cairocffi as cairo

""" Set of classes for helping deal with SVG graphics """


class SVG():
    """
    Class is a container for basic drawing objects.  It is used to record
    actions that make up an svg so that they may be redrawn using Jquery, pycairo
    or another drawing tool
    """

    def __init__(self):
        self.operations = []

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

    def render_png(self, writer, width, height):
        """
        Renders this svg object as a PNG
        @param writer - any file like object that has a write method
        @param width - width of the png
        @param height - height of the png
        """

        surface = cairo.ImageSurface (cairo.FORMAT_ARGB32, width, height)
        context = cairo.Context (surface)

        for op in self.operations:
            if op.type == 'line':
                context.move_to(op.x+.5, op.y+.5)
                context.line_to(op.x1+.5, op.y1+.5)
                context.set_line_width(op.stroke)
                r,g,b = RGBTuple(op.color)
                context.set_source_rgba(r,g,b,1)
                context.stroke()

            elif op.type == 'text':
                red, green, blue = RGBTuple(op.color)
                context.set_source_rgba(red,green,blue,1)
                context.set_font_size (op.size)
                context.move_to (op.x, op.y)
                context.show_text (op.text)

            elif op.type == 'rect':
                context.rectangle(op.x, op.y, op.width, op.height)
                if op.color and op.color <> 'None':
                    red, green, blue = RGBTuple(op.color)
                    context.set_source_rgba(red,green,blue,1)
                    context.set_line_width(op.stroke)
                    context.stroke_preserve()
                if op.fill and op.fill <> 'None':
                    r,g,b = RGBTuple(op.fill)
                    context.set_source_rgba(r,g,b,1)
                else:
                    context.set_source_rgba(0,0,0,0)
                context.fill()

        surface.write_to_png(writer)



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


def RGBTuple(rgbString):
    """
    Converts a hex string to a tuple of RGB integer values
    """
    sub = rgbString[-6:]
    red = int(sub[:2],16)/255.0
    green = int(sub[2:4],16)/255.0
    blue = int(sub[4:], 16)/255.0
    return (red,green,blue)