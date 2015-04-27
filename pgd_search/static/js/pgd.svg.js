function render_svg(svg, paper, font, func) {
    /*
     renders a pgd.svg object.  PGD uses a custom SVG like object that contains
     a list of drawing operations.  It is converted to javascript lists and
     dictionaries using json.  This function iterates through the list of
     operations and processes them.
     
     Rectangles can supply additional data with the objects.  This allows
     objects to have user interactions such as mouseover, click, etc.  A user
     supplied function will process the additional data and setup whatever
     additional logic is required
     
     @param svg - dictionary of svg drawing operations
     @param paper - Raphael paper (canvas) object
     @param font - Raphael font object
     @param func - function used to process additional data that is
                   included with rectangles.
    */
    
    var aafix = 0.5;
    var rects_width_fix = 0;
    var rects_height_fix = 0;
    var line_x_aafix = 0.5;
    var line_y_aafix = 0.5;
    var rects_x_aafix = 0;
    var rects_y_aafix = 0;
    var is_IE = false;
    
    // IE checks - check different things not supported by IE to detect them
    if($.jqbrowser.msie()){
        is_IE = true;
        rects_width_fix = 0;
        rects_height_fix = 0;
        line_x_aafix = 2;
        line_y_aafix = 2;
        rects_x_aafix = 0.5;
        rects_y_aafix = 0.5;
    } else if ($.jqbrowser.win()) {
        aafix = 0;
        rects_x_aafix = -1;
        rects_y_aafix = -.5;
        line_x_aafix = -.5;
        line_y_aafix = 0;
        rects_width_fix = 1.5;
    }
    
    for (i=0; i<svg.length; i++) {
        op = svg[i]
        if (op['type'] == 'line') {
            l = paper.path("M{0},{1}L{2},{3}",
                           op['x']+line_x_aafix,
                           op['y']+line_x_aafix,
                           op['x1']+line_y_aafix,
                           op['y1']+line_y_aafix)
                .attr({stroke: op['color']});
        } else if (op['type'] == 'text') {
            t = paper.print(op['x'],op['y'],op['text'], font, op['size']);
            if (op['rotate'] != 0) {
                t.rotate(op['rotate'], op['x'],op['y']);
            }
        } else if (op['type'] == 'rect') {
            r = paper.rect(op['x']+rects_x_aafix,
            op['y']+rects_y_aafix,
            op['width']+rects_width_fix,
            op['height']+rects_height_fix);
            if(op['fill'] != undefined) {
                r.attr('fill', op['fill']);
            }
            r.attr('stroke', op['color'])
            r.attr('stroke-width',op['stroke']);
            
            if (op['data'] != undefined) {
                func(r, op);
            }
        }
    }
}