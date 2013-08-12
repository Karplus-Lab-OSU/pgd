n = Math.ceil(MAX_LENGTH / 2);
iArray = new Array(MAX_LENGTH - 1);
for (i=0; i<n; i++) {
    iArray[i*2] = i+1;
    if (i*2+1 < MAX_LENGTH - 1){
        iArray[i*2+1] = 0-i-1;
    }
}

function getNegateRow(target) {
    return negateRow[$(target).parent().attr('class').substring(10,12)];
}

function putNegateRow(target, value) {
    negateRow[$(target).parent().attr('class').substring(10,12)] =
        value;
}

function sizeChange(evt){
    oldValue = parseInt($('#oldLength').val());
    newValue = parseInt($('#id_residues').val());
    if (oldValue != newValue) {
        if(qtipTotal!=0){
            qtipTotal = 0;
            $('div.qtip').qtip("destroy");
            createQtips();
        }

        //removing values
        if (oldValue > newValue) {
            for (i = MAX_LENGTH - 1; i >= newValue; i--) {
                str = '.col_' + iArray[i-1]

                // TODO: I think this code might be getting run too many times over.
                fieldsToReset = document.getElementsByClassName("col_"+iArray[i-1]);
                console.log(fieldsToReset.length+"  fields to reset");
                for (var i_col=0; i_col<fieldsToReset.length; i_col++) {
                    // Hack for #9495: Uncheck all selections in removed columns.
                    inputArray = fieldsToReset[i_col].getElementsByClassName("needs_reset");
                    for (var k=0; k<inputArray.length; k++) {
                        if (inputArray[k].type == "checkbox") {
                            inputArray[k].checked = false;   // Uncheck checkbox.
                            $(inputArray[k]).parent().attr("class", ' ');   // Classname is set to an empty space here because that's what happens when residue is unchecked by hand. This artifact should be figured out, but I don't want to right now. TODO.
                        }
                        else {
                            inputArray[k].value = '';   // Eh.. this is so that isDefault() inside updateInclusionField() will see that the field style of this element should be reset.
                            updateInclusionField(inputArray[k], false);   // Reset field style.
                        }
                    }
                }

                $(str).hide();
                $(str + ' input.field').attr("disabled", true);
            }
        // adding values
        } else {
            for(i=0; i<=newValue; i++){
                str = '.col_' + iArray[i-2]
                $(str + ' input.field').attr("disabled", false);
                $(str).show();
                /* Hax for #1565, #1566, and #9495: Add the
                 * initial angles here instead of on the form.
                 * This is because the forms are super-wimpy and
                 * don't store state very well. */
                if (evt != 1) {   // Hax (maybe?) for #9495. sizeChange() is called when the document is readied with an argument of 1.
                    $("#id_ome_" + iArray[i - 2]).val("<=-90,>=90");
                    frobInclusions($("#id_ome_" + iArray[i - 2]));
                    $("#id_bm_" + iArray[i - 2]).val("<25");
                    frobInclusions($("#id_bm_" + iArray[i - 2]));
                    $("#id_bg_" + iArray[i - 2]).val("<25");
                    frobInclusions($("#id_bg_" + iArray[i - 2]));
                    $("#id_bs_" + iArray[i - 2]).val("<25");
                    frobInclusions($("#id_bs_" + iArray[i - 2]));
                }
            }
        }
        //update length
        $('#oldLength').val(newValue);

        //update size of table
        //if (newValue < 5){
        //    val = 54;
        //} else {
            val = 10.3*newValue+6
        //}

        // Update things to include in search. Fixes #9495.
        $('tr.peptide_row .multiselect li, tr.conf_row .multiselect li').parent().updateInclusionSelect();

        $('table#residues').css({'width':val+'em'});
        $('table#protein').css({'width':600});
    }
}

function validateField(val) {
    /* Clean up field by removing whitespace and replacing OR
     * synonyms with commas. */
    var regex = /\s+/g;
    var val = val.replace(regex, '');
    regex = /(\|){2}|(or){1}/g;
    val = val.replace(regex, ',');

    if (val == '') {
        return true;
    }

    var m = match(val);
    // Apparently empty lists are returned as empty strings
    if (m != '') {
        return true;
    }
    return false;
}

/* Update the classes on and around any node which can be an
 * inclusion input of some sort, indicating if the node is an
 * inclusion or exclusion filter. This works on both fields and
 * selections. */
function frobInclusions(node) {
    var $node = $(node);
    var $parent = $node.parent();
    var val = $parent.children(".include").val();
    if (val == '') {
        val = 1;
    }
    if (val == 1 && !$node.hasClass('field_include')) {
        $parent.children(".toggle").addClass('toggle_Include');
        $parent.children(".include").val(1);
        $parent.children(".include").attr('disabled', false);
        $node.addClass('field_include');
        $node.removeClass('field_exclude');
        putNegateRow($parent, true);
    } else if (val == 0 && !$node.hasClass('field_exclude')) {
        $parent.children(".toggle").addClass('toggle_Exclude');
        $parent.children(".include").val(0);
        $parent.children(".include").attr('disabled', false);
        $node.addClass('field_exclude');
        $node.removeClass('toggle_Include');
        putNegateRow($parent, false);
    }
}

function updateInclusionField(node, syntax) {
    var $node = $(node);
    var $parent = $node.parent();
    var val = $node.val();

    // if field is not null check for syntax errors
    if (syntax) {
        if (validateField(val)) {
            $node.removeClass('field_error');
        } else {
            $parent.children(".toggle").removeClass('toggle_Exclude toggle_Include');
            $node.addClass('field_error');
            $node.removeClass('field_exclude field_include');

            /* Return immediately so only the error class is
             * applied. */
            return;
        }
    }

    // set field style based on if the field has the default value
    if (isDefault(node, syntax)) {
        $parent.children(".toggle").removeClass('toggle_Include toggle_Exclude');
        $parent.children(".include").val('');
        $parent.children(".include").attr('disabled', true);
        $node.removeClass('field_exclude field_include');
    } else {
        frobInclusions(node);
    }
};

function updateInclusionSelect(node) {
    var $node = $(node);
    var $parent = $node.parent();
    var is_checked = false;

    // Is anything checked?
    $node.find('input').each(function() {
        if (this.checked) {
            is_checked = true;
        }
    });

    if (is_checked) {
        frobInclusions(node);
    } else {
        $node.removeClass('field_exclude field_include');
        $parent.children(".toggle").removeClass('toggle_Include toggle_Exclude');
        $parent.children(".include").val('');
    }
};

// extend function to jquery
$.fn.updateInclusionField = function(){
    this.each(function(){
        updateInclusionField(this, true);
    });
};
$.fn.updateInclusionSelect = function(){
    this.each(function(){
        updateInclusionSelect(this);
    });
};

/*
Determine if an element has the default value for that field type
TODO: The name of this function is misleading, because it currently just checks if the field is blank (disregarding whitespace). A better name might be something like "isEmpty".
*/
function isDefault(element, field){
    val = $(element).val();
    // strip whitespace for input boxes
    if (field) {
        val = val.replace(/\s+/g, '');
        return val == '';
    }

    // select boxes just check for null
    return val == '' || val == null;
};

function createQtips(){
    $('#residuesContainer input.field,#residuesContainer td.field > .multiselect')
    .not('.sidechain_length_input input,.sidechain_angle_input input')//Don't make qtips for sidechain angle input boxes
    .each(function(){
        currentID = $(this).attr('id');
        col = currentID.substring(currentID.length-1,currentID.length);
        sign = currentID.substring(currentID.length-2,currentID.length-1);
        if(col==Math.ceil(parseInt($('#id_residues').val()/2)) && col!= null && col!='' && sign!='-'){
            qtipTotal=qtipTotal+1;
            currentID = currentID.substring(0, currentID.length - 1 );
            target = currentID.replace('-','');
            target = ('#' + target + col);
            helpName = currentID.substring(3, currentID.length - 1 );
            helpName = helpName.replace('_','');
            helpName = helpName.replace('_','');
            helpName = helpName.replace('_','');
            helpName = '#' + helpName.replace(' ','');
            contentText = $(helpName).html();
            titleText = $(helpName).attr('name');
            $(this).qtip({
            position:{
                adjust: { //Make the tooltip avoid going outside the window due to browser resizing
                        resize: true
                        },
                corner: {
                        target: 'rightMiddle',
                        tooltip: 'leftMiddle'
                        },
                target: $(target)
                     },
            style:{ //Style the qtip to look like the page
                  color: 'black',
                  textAlign: 'center',
                  title: { color: 'white', background: '#666666'},
                  width: { min: 350, max: 430 },
                  border:{ 
                         width: 6, 
                         color: '#909090'
                         },
                  tip: { //Options for the tail
                       corner: 'leftMiddle',
                       color: '#909090',
                       size: { x: 20, y : 15 }
                       }
                   },
            show: false,
            hide: false,
            content:{
                    prerender: false,//This kills pageload times if enabled
                    text: contentText,
                    title:  {
                            text: titleText,
                            button: 'Close' //Enables a button that closes the qtip
                            }
                    }
                });
            }
        });
    $('#protein td>select,#protein td>input').not('#id_proteins_i,#oldLength').each(function(){
            currentID = $(this).attr('id');
            currentID = currentID.replace(' ','');
            if(currentID!='' && currentID!=null && currentID!=undefined){
                if(currentID.substring(currentID.length-3,currentID.length)=='Min'){
                    target='#'+currentID.replace('Min','Max');
                }
                else{
                    target='#'+currentID;
                }
                helpName = '#'+currentID.replace('_','');
                contentText = $(helpName).html();
                titleText = $(helpName).attr('name');
            }
            qtipTotal=qtipTotal+1;
            $(this).qtip({
                position:{
                            adjust: { //Make the tooltip avoid going outside the window due to browser resizing
                                    resize: true
                                    },
                            corner: {
                                    target: 'rightMiddle',
                                    tooltip: 'leftMiddle'
                                    },
                            target: $(target)
                         },
                style:{ //Style the qtip to look like the page
                      color: 'black',
                      textAlign: 'center',
                      title: { color: 'white', background: '#666666'},
                      width: { min: 350, max: 430 },
                      border:{
                             width: 6,
                             color: '#909090'
                             },
                      tip: { //Options for the tail
                           corner: 'leftMiddle',
                           color: '#909090',
                           size: { x: 20, y : 15 }
                           }
                       },
                show: false,
                hide: false,
                content:{
                        prerender: false,//This kills pageload times if enabled
                        text: contentText,
                        title:  {
                                text: titleText,
                                button: 'Close' //Enables a button that closes the qtip
                                }
                        }
                });
        });
};


function create_angles_helper( qcontent, qtarget, anch ){
    $(anch).qtip({
        position:{
            adjust: { //Make the tooltip avoid going outside the window due to browser resizing
                    resize: true
                    },
            corner: {
                    target: 'topMiddle',
                    tooltip: 'bottomMiddle'
                    },
            target: qtarget
                 },
            style:{ //Style the qtip to look like the page
                color: 'black',
                textAlign: 'center',
                title: { color: 'white', background: '#666666'},
                width: { min: 200, max: 200 },
                border:{ 
                        width: 6, 
                        color: '#909090'
                        },
                tip: { //Options for the tail
                    corner: 'bottomMiddle',
                    color: '#909090',
                    size: { x: 20, y : 15 }
                    },
                classes:{
                    content: "angle_visualizer",
                    tooltip: "angle_window"
                    }
                },
            show:false,
            hide:false,
            content:{
                    prerender: false,//This kills pageload times if enabled
                    text: qcontent,
                    title:  {
                            text: "Angle Visualizer",
                            }
                    }
            })
};

function toggleQtips(){
    if(qtipsOn){
        $('div.qtip').qtip('hide');
        $('div.qtip').qtip("disable");
        qtipsOn=false;
        $('.qtipToggle').attr('value','Turn Help On');
    }
    else{
        $('div.qtip').qtip("enable");
        qtipsOn=true;
        $('.qtipToggle').attr('value','Turn Help Off');
    }

    qtipState = '';
};

var r = undefined;
var _angles = undefined;
var rad = Math.PI/180;
var cur = undefined;

/* Draw an arc on the canvas.
 * cx is the X coordinate of the circle containing this arc.
 * cy is the Y coordinate of the circle containing this arc.
 * start is the beginning angle of the arc.
 * end is the ending angle of the arc.
 * R is the radius of the arc.
 */
function draw_arc(cx, cy, start, end, R) {
    var x1 = cx - R * Math.cos(-start * rad),
    x2 = cx - R * Math.cos(-end * rad),
    y1 = cy - R * Math.sin(-start * rad),
    y2 = cy - R * Math.sin(-end * rad);
    if (-180==start && 180==end) {
        path = ["M", 70.001, 90.001 - R, "A", R, R, 0, 1, 1, 69.99, 90.001 - R].join();
    } else {
        path = ["M",cx,cy,"L", x1, y1, "A", R, R, 0, +(end - start > 180), 0, x2, y2, "z"].join();
    }
    return r.path(path).attr({stroke: '#004400', fill:'#c2e192', "stroke-width": 1, opacity:.8});
}

/* Generate arcs for lists of angles. */
function angles(list) {
    if (_angles != undefined) {
        for (i in _angles) {
            _angles[i].remove();
        }
    }
    base = 100;
    _angles = [];
    /* If the list of angles is empty, push a circle which
     * encompasses the entire range. */
    if (list.length == 0) {
        _angles.push(r.circle(80, 90, 40).attr({stroke: '#004400', fill:'#c2e192', "stroke-width": 1, opacity:.8}));
    } else {
        for (i in list) {
            _angles.push(draw_arc(80,90,list[i][0]+180,list[i][1]+180, 40));
        }
    }
}

//Checks the list of numbers for compliance to search ranges
function check_range(val) {
    var vals = null;

    for (range in val) {
        vals = val[range];
        // Must go counterclockwise around the circle...
        if ((vals[0] - vals[1]) >= 0) {
            return false;
        }
        //Must be in the range of -180 to 180
        if (vals[0]> 180 || vals[0] < -180) {
            return false;
        }
        if (vals[1]> 180 || vals[1] < -180) {
            return false;
        }
    }
    return true;
}

/* Attempt to parse a float with parseFloat, raising an exception
 * if the string couldn't be parsed. */
function parseFloatOrFail(s) {
    var retval = parseFloat(s);
    if (isNaN(retval)) {
        throw "Couldn't parse float" + s;
    }
    return retval;
}

function match(val) {
    var vals = val.split(',');
    var processed = [];

    try {
        for (i in vals) {
            v = vals[i];
            if (v.indexOf('<=') != -1) {
                processed.push([-180, parseFloatOrFail(v.substring(2))]);
            } else if (v.indexOf('>=') != -1) {
                processed.push([parseFloatOrFail(v.substring(2)), 180]);
            } else if (v.indexOf('<') != -1) {
                processed.push([-180, parseFloatOrFail(v.substring(1))]);
            } else if (v.indexOf('>') != -1) {
                processed.push([parseFloatOrFail(v.substring(1)), 180]);
            } else if (v.indexOf('-',1) != -1) {
                z = v.indexOf('-', 1);
                v1 = parseFloatOrFail(v.substring(0, z));
                v2 = parseFloatOrFail(v.substring(z + 1));
                processed.push([v1, v2]);
            } else {
                v = parseFloatOrFail(v);
                processed.push([v, v + 1]);
            }
        }
    } catch (error) {
        return [];
    }

    if (check_range(processed)) {
        return processed;
    } else {
        return [];
    }
}

function angles_helper_init(){
    r = Raphael("canvas", 170, 175);
    param = {stroke: "#000", "stroke-width": 1, "stroke-dasharray":"- "}
    
    x = 80;
    y = 90;
    radius = 30;
    r.circle(x,y,radius).attr({stroke: "#AAAAAA", "stroke-width": 1, "stroke-dasharray":"- "});
    
    r.path("M"+x+" "+(y-radius-40)+"L"+x+ " " +(y+radius+40)).attr(param);
    r.path("M"+(x-radius-50)+" "+y+"L"+(x+radius+40)+ " " +y).attr(param);
    r.text(x+radius+48, y, "0");
    r.text(x, y-radius-48, "90");
    r.text(x, y+radius+45, "-90");
    r.text(x-radius-35, y-10, "180");
    r.text(x-radius-35, y+12, "-180");
};

function redraw(){
    
    //Add or remove sidechain angles
    active_sidechains = [];
    for (active_aa in selected_aa){//Add all of the sidechains into a list using the active aa's and the lookup table
        sidechains = sidechain_angle_lookup[selected_aa[active_aa]];
        for (sidechain in sidechains){
            active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
        }
    }
    header = document.getElementById('sidechain_angle_header');
    if (header.innerHTML == 'Sidechain Angles ↓'){//Only show if header is open.
        for (angle in active_sidechains){
            var row = $("tr#" + active_sidechains[angle] + "_row");
            row.show();
            row.children('td').children('span').children('input').hide();
            if (protList[active_sidechains[angle].substring(0,3)]){
                row.children('td').each(
                    function(){
                        //CHANGE HERE FOR NEGATION CLAUSE
                        if (getNegateRow(this) != false && protList[active_sidechains[angle].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                            //re-enable the angles if they apply to the column.
                            $(this).children('span').children('input').show();
                        }
                    }
                    
                );
            }
        }
    }
    
    for (i in AAs){//Only remove those sidechains that aren't still shown by another aa
        aa = AAs[i];
        if (selected_aa.indexOf(aa) == -1) {
            for (angle in sidechain_angle_lookup[aa]){
                $('tr#'+aa+'__'+sidechain_angle_lookup[aa][angle]+'_row').hide();
            }
            for (sclength in sidechain_length_lookup[aa]){
                $('tr#'+aa+'__'+sidechain_length_lookup[aa][sclength]+'_row').hide();
            }
        }
    }


    /*********************************************************/
    //Add or remove sidechain lengths
    active_sidechains = [];
    for (active_aa in selected_aa){//Add all of the sidechains into a list using the active aa's and the lookup table
        sidechains = sidechain_length_lookup[selected_aa[active_aa]];
        for (sidechain in sidechains){
            active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
        }
    }
    //unfiltered_sidechains = sidechain_length_lookup[aa];
    header = document.getElementById('sidechain_length_header');
    if (header.innerHTML == 'Sidechain Lengths ↓'){//Only show if header is open.
        for (sclength in active_sidechains){
            var row = $("tr#" + active_sidechains[sclength] + "_row");
            row.show();
            row.children('td').children('span').children('input').hide();
            if (protList[active_sidechains[sclength].substring(0,3)]){
                row.children('td').each(
                    function(){
                        //access to form element via $(this)
                        
                        if (getNegateRow(this) != false && protList[active_sidechains[sclength].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                            //re-enable the angles if they apply to the column.
                            $(this).children('span').children('input').show();
                        }
                    }
                );
            }
        }
    }
    
};

function fullRedraw(aa, aaColumn){
    /*********************************************************/
    //Add or remove sidechain angles
    active_sidechains = [];
    for (active_aa in selected_aa){//Add all of the sidechains into a list using the active aa's and the lookup table
        sidechains = sidechain_angle_lookup[selected_aa[active_aa]];
        for (sidechain in sidechains){
            active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
        }
    }
    header = document.getElementById('sidechain_angle_header');
    if (header.innerHTML == 'Sidechain Angles ↓'){//Only show if header is open.
        for (angle in active_sidechains){
            var row = $("tr#" + active_sidechains[angle] + "_row");
            row.show();
            row.children('td').children('span').children('input').hide();

            if (protList[active_sidechains[angle].substring(0,3)]){
                row.children('td').each(
                    function(){
                        //access to form element via $(this)
                        if (getNegateRow(this) != false && protList[active_sidechains[angle].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                            //re-enable the angles if they apply to the column.
                            $(this).children('span').children('input').show();
                        }
                    }
                    
                );
            }
        }
    }

    // Clear and remove all sidechain angles not shown by another aa.
    if((!$this.hasClass('selected')) && ((selected_aa.indexOf(aa)==-1))){
        for (angle in sidechain_angle_lookup[aa]){
            fieldToReset = $('input#id_'+aa+'__'+sidechain_angle_lookup[aa][angle]+'_'+aaColumn);
            fieldToReset.val('');   // Clear field value.
            updateInclusionField(fieldToReset, false);   // Update field style.

            $('tr#'+aa+'__'+sidechain_angle_lookup[aa][angle]+'_row').hide();
        }
    }
    /*********************************************************/
    //Add or remove sidechain lengths
    active_sidechains = [];
    for (active_aa in selected_aa){//Add all of the sidechains into a list using the active aa's and the lookup table
        sidechains = sidechain_length_lookup[selected_aa[active_aa]];
        for (sidechain in sidechains){
            active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
        }
    }
    //unfiltered_sidechains = sidechain_length_lookup[aa];
    header = document.getElementById('sidechain_length_header');
    if (header.innerHTML == 'Sidechain Lengths ↓'){//Only show if header is open.
        for (sclength in active_sidechains){
            var row = $("tr#" + active_sidechains[sclength] + "_row");
            row.show();
            row.children('td').children('span').children('input').hide();
            if (protList[active_sidechains[sclength].substring(0,3)]){
                row.children('td').each(
                    function(){
                        //access to form element via $(this)
                        if (getNegateRow(this) != false && protList[active_sidechains[sclength].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                            //re-enable the angles if they apply to the column.
                            $(this).children('span').children('input').show();
                        }
                    }
                    
                );
            }
        }
    }

    // Clear and remove all sidechain lengths not shown by another aa.
    if((!$this.hasClass('selected')) && ((selected_aa.indexOf(aa)==-1))){
        for (sclength in sidechain_length_lookup[aa]){
            fieldToReset = $('input#id_'+aa+'__'+sidechain_length_lookup[aa][sclength]+'_'+aaColumn);
            fieldToReset.val('');   // Clear field value.
            updateInclusionField(fieldToReset, false);   // Update field style.

            $('tr#'+aa+'__'+sidechain_length_lookup[aa][sclength]+'_row').hide();
        }
    }
    /*********************************************************/
};

$(document).ready(function() {
    $('#helpPage').hide();
    
    if(qtipsOn){
        createQtips();
    }
    angles_helper_init();
    $('#canvas').hide();//Hide this after the init, qtip will grab it once the qtip is made
    $('span.toggle').click(function(evt){
        var $tar = $(evt.target);
        if ($tar.hasClass('toggle_Include')) {
            $tar.removeClass('toggle_Include');
            $tar.addClass('toggle_Exclude');
            putNegateRow($tar.parent(), false);
            redraw();
            $tar.parent().children(".include").val(0);
            $tar.parent().children(".field, #id_proteins").addClass('field_exclude').removeClass('field_include');
        }else if ($tar.hasClass('toggle_Exclude')){
            $tar.removeClass('toggle_Exclude');
            $tar.addClass('toggle_Include');
            putNegateRow($tar.parent(), true);
            redraw();
            $tar.parent().children(".include").val(1);
            $tar.parent().children(".field, #id_proteins").addClass('field_include').removeClass('field_exclude');
        }
    });

        $('input.field')
        .keyup(function(evt){
            updateInclusionField(evt.target, true);
            //update svg
            val = $(this).val();
            m = match(val);
            angles(m);
            })
        .change(function(evt){
                updateInclusionField(evt.target, true);
            })
        .not('.lengths_row input, mobility_row input, .sidechain_length_input input')
            .blur(function(evt){
                angle_qtip_api.hide();
            })
            .focus(function(evt){
                //update svg
                val = $(this).val();
                m = match(val);
                angles(m);
                if(firstRun == true){//initial run needs to have target set and qtip made
                    create_angles_helper('Text',$('#'+evt.target.id),'#'+evt.target.id)
                    angle_qtip_api = $("#"+evt.target.id).qtip('api');
                    angle_qtip_api.show();
                    $('#canvas').show();//show it so the qtip can grab it
                    $('.angle_visualizer').html($('#canvas'));
                    firstRun = false;
                }
                angle_qtip_api.options.position.target = $('#'+evt.target.id);
                //Set a slight delay to prevent the hide function from stepping on the show function
                window.setTimeout(function(){angle_qtip_api.updatePosition($('#'+evt.target.id),false).show();}, 200);
            }).updateInclusionField();


    $('#protein #id_proteins')
        .autocomplete(SITE_ROOT+'/search/protein_code/', {'multiple':true, 'autoFill':true})
        .keyup(function(evt){
                updateInclusionField(evt.target, false);
            })
        .change(function(evt){
                updateInclusionField(evt.target, false);
            })
        .each(function(){updateInclusionField(this, false);});
    $('#id_residues').change(sizeChange);

    //force update of segments
    $('#oldLength').val(1);
    sizeChange(1);   // TODO: What does the 1 do here?

    // setup multiselects
    $('tr.peptide_row .multiselect li, tr.conf_row .multiselect li')
        .click(function(){
            $this = $(this);
            aa = $this.html().substring(0,3).toUpperCase();
            aaColumn = $this.children('input').attr('name').substring(3,5);
            if ($this.hasClass('selected')){
                $this.removeClass('selected');
                $this.children('input').attr('checked',false);
                if (protList[aa] != undefined){//If a Dict entry not exists, add new column
                    protList[aa][aaColumn] = false;
                    
                }
                else{//Otherwise add it to Dict and add new column
                    protList[aa] = [];
                    protList[aa][aaColumn] = false;
                    
                }
            }else{
                $this.addClass('selected');
                $this.children('input').attr('checked',true);
                if (protList[aa] != undefined){//If a Dict entry not exists, add new column
                    protList[aa][aaColumn] = true;
                    
                }
            }
            $this.parent().updateInclusionSelect()
        })
        // update styles based on if its already checked
        .each(function(){
            $this = $(this);
            if ( $this.children('input').attr('checked')){
                 $this.addClass('selected');
            }
        });

    
    $('tr.peptide_row .multiselect li.selected')//Add any selected AAs left over from a previous state
        .each(function(){
            $this = $(this);
            aa = $this.html().substring(0,3).toUpperCase();
            aaColumn = $this.children('input').attr('name').substring(3,5);
            selected_aa.push(aa);
            if (protList[aa] != undefined){//If a Dict entry not exists, add new column
                protList[aa][aaColumn] = true;
                
            }
            else{//Otherwise add it to Dict and add new column
                protList[aa] = [];
                protList[aa][aaColumn] = true;
                
            }
            
        });

    // setup aa multiselects to work with sidechain show/hide

    $('tr.peptide_row .multiselect li')
        .click(function(){
            $this = $(this);
            //Strip out aa and manage list of active aa's
            if ($this.html() != null){
                aa = $this.html().substring(0,3).toUpperCase();
            }
            
            if ($this.children('input').attr('name') != null){
                aaColumn = $this.children('input').attr('name').substring(3,5);
            }
            
            if ($this.hasClass('selected')){
                selected_aa.push(aa);
                
                if (protList[aa] != undefined){//If a Dict entry not exists, add new column
                    protList[aa][aaColumn] = true;
                    
                }
                else{//Otherwise add it to Dict and add new column
                    protList[aa] = [];
                    protList[aa][aaColumn] = true;
                    
                }
            }
            else{
                aa_to_remove = selected_aa.indexOf(aa)
                selected_aa.splice(aa_to_remove,1);
                
                if (protList[aa]){//If a Dict entry not exists, add new column
                    protList[aa][aaColumn] = false;
                }
            }
            fullRedraw(aa, aaColumn);
        
        });
    

    $('tr.peptide_row .multiselect, tr.conf_row .multiselect').updateInclusionSelect();

    $('tr.chi_row').hide();
    $('#chi_header').click(function(evt){
        var header = document.getElementById('chi_header');
        if (header.innerHTML == 'χ Angles →'){//Must be χ, it won't match &chi;
            header.innerHTML = 'χ Angles ↓';
            $('tr.chi_row').show();
        } else {
            header.innerHTML = 'χ Angles →';
            $('tr.chi_row').hide();
        }
    });

    $('tr.sidechain_angle_row').hide();
    $('#sidechain_angle_header').click(function(evt){
        var header = document.getElementById('sidechain_angle_header');
        if (header.innerHTML == 'Sidechain Angles →'){
            header.innerHTML = 'Sidechain Angles ↓';
            active_sidechains = [];
            for (active_aa in selected_aa){
                sidechains = sidechain_angle_lookup[selected_aa[active_aa]];
                for (sidechain in sidechains){
                    active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
                }
            }
            for (sidechain in active_sidechains){
                $('tr#'+active_sidechains[sidechain]+'_row').show();
                //set all cells to disabled as default; cells are more likely to need disabling than enabling
                $('tr#'+active_sidechains[sidechain]+'_row').children('td').children('span').children('input').hide();
                if (protList[active_sidechains[sidechain].substring(0,3)]){
                    
                    $('tr#'+active_sidechains[sidechain]+'_row').children('td').each(
                        function(){
                            //access to form element via $(this)
                            if (getNegateRow(this) != false && protList[active_sidechains[sidechain].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                                //re-enable the angles if they apply to the column.
                                $(this).children('span').children('input').show();
                            }
                        }
                        
                    );
                }
                //alert(active_sidechains[sidechain]);
            }
        } else {
            header.innerHTML = 'Sidechain Angles →';
            $('tr.sidechain_angle_row').hide();
            
        }
    });

    $('tr.sidechain_length_row').hide();
    $('#sidechain_length_header').click(function(evt){
        var header = document.getElementById('sidechain_length_header');
        if (header.innerHTML == 'Sidechain Lengths →'){
            header.innerHTML = 'Sidechain Lengths ↓';
            active_sidechains = [];
            for (active_aa in selected_aa){
                sidechains = sidechain_length_lookup[selected_aa[active_aa]];
                for (sidechain in sidechains){
                    active_sidechains.push(selected_aa[active_aa]+'__'+sidechains[sidechain]);
                }
            }
            for (sidechain in active_sidechains){
                $('tr#'+active_sidechains[sidechain]+'_row').show();
                //set all cells to disabled as default; cells are more likely to need disabling than enabling
                $('tr#'+active_sidechains[sidechain]+'_row').children('td').children('span').children('input').hide();
                if (protList[active_sidechains[sidechain].substring(0,3)]){
                    
                    $('tr#'+active_sidechains[sidechain]+'_row').children('td').each(
                        function(){
                            //access to form element via $(this)
                            
                            if (getNegateRow(this) != false && protList[active_sidechains[sidechain].substring(0,3)][$(this).attr('class').substring(10,12)] == true){
                                //re-enable the angles if they apply to the column.
                                $(this).children('span').children('input').show();
                            }
                        }
                        
                    );
                }
            }
        } else {
            header.innerHTML = 'Sidechain Lengths →';
            $('tr.sidechain_length_row').hide()
            
        }
    });

    $('#peptide_header').click(function(evt){
        var header = document.getElementById('peptide_header');
        if (header.innerHTML == 'Composition →'){
            header.innerHTML = 'Composition ↓';
            $('tr.peptide_row').show();
        } else {
            header.innerHTML = 'Composition →';
            $('tr.peptide_row').hide();
        }
    });
    
    $('#conf_header').click(function(evt){
        var header = document.getElementById('conf_header');
        if (header.innerHTML == 'Conformation →'){
            header.innerHTML = 'Conformation ↓';
            $('tr.conf_row').show();
        } else {
            header.innerHTML = 'Conformation →';
            $('tr.conf_row').hide();
        }
    });

    $('tr.mobility_row').hide();
    $('#mobility_header').click(function(evt){
        var header = document.getElementById('mobility_header');
        if (header.innerHTML == 'Mobility →'){
            header.innerHTML = 'Mobility ↓';
            $('tr.mobility_row').show();
        } else {
            header.innerHTML = 'Mobility →';
            $('tr.mobility_row').hide();
        }
    });

    $('tr.angles_row').hide();
    $('#angles_header').click(function(evt){
        var header = document.getElementById('angles_header');
        if (header.innerHTML == 'Angles →'){
            header.innerHTML = 'Angles ↓';
            $('tr.angles_row').show();
        } else {
            header.innerHTML = 'Angles →';
            $('tr.angles_row').hide();
        }
    });

    $('tr.lengths_row').hide();
    $('#lengths_header').click(function(evt){
        var header = document.getElementById('lengths_header');
        if (header.innerHTML == 'Lengths →'){
            header.innerHTML = 'Lengths ↓';
            $('tr.lengths_row').show();
        }
        else{
            header.innerHTML = 'Lengths →';
            $('tr.lengths_row').hide();
        }
    });

    $('#residuesContainer input.field,#residuesContainer td.field>.multiselect').not('.sidechain_length_input input,.sidechain_angle_input input').mousedown(function(event){
            $clickedOn  = $(event.target);
            if($clickedOn.is('input')){
                rowTarget = event.target.id;
            }
            else if($clickedOn.is('li')){
                rowTarget = event.target.parentNode.id;
            }
            if(rowTarget!=''&&rowTarget!=undefined){
            rowTarget = rowTarget.substring(0, rowTarget.length-1);
            rowTarget = rowTarget.replace('-','');
            rowTarget = '#' + rowTarget + Math.ceil(parseInt($('#id_residues').val()/2));
            rowTarget = rowTarget.replace(' ','');
            if(rowTarget!=qtipState&&qtipsOn){
                $('div.qtip').not('div.angle_window').qtip('hide');//XXX Exclude angle visualizer to avoid flickering
                $(rowTarget).not('div.angle_window').qtip('show');
            }
            qtipState=rowTarget;
            rowTarget='';
            }
        });

    $('#protein td>select,#protein td>input').not('#id_proteins_i,#oldLength').mousedown(function(event){
        $clickedOn  = $(event.target);
        if($clickedOn.is('input')){
            rowTarget =event.target.id;
        }
        else if($clickedOn.is('select')){
            rowTarget =event.target.id;
        }
        if(rowTarget!=''&&rowTarget!=undefined&&qtipsOn&&rowTarget!=qtipState){
            $('div.qtip').not('div.angle_window').qtip('hide');
            $('#'+rowTarget).not('div.angle_window').qtip('show');
            qtipState=rowTarget;
            rowTarget='';
        }
    });
    
    $('#search').submit(function(){
        /* Disable all unused parameters so the data is not sent */
        $this = $(this);
        $this.find('input').each(function(){
            if ($(this).val() == '') {
                $(this).attr('disabled',true);
            }
        });
    });
});

// vim: expandtab

