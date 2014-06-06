from django.template.defaultfilters import stringfilter, floatformat
from django.utils.safestring import mark_safe
from django import template
register = template.Library()

from pgd_constants import AA_CHOICES_DICT
from pgd_search.plot.PlotForm import PROPERTY_CHOICES_DICT


ROMAN_TO_GREEK = dict(
    A='&alpha;',
    B='&beta;',
    G='&gamma;',
    D='&delta;',
    E='&epsilon;',
    H='&eta;',
    Z='&zeta;'
)


@register.filter(name='full_aa')
@stringfilter
def full_aa(value):
    """
    Filter that returns the full AA code from the 1 letter AA code
    """
    return AA_CHOICES_DICT[value]
register.filter('full_aa',full_aa)


@register.filter(name='index_lookup')
def index_lookup(value, arg):
    """
    Filter that allows you to lookup the value of a list using another variable in the page
    """
    try:
        return value[int(arg)] if isinstance(value, (list,tuple)) else value[arg]
    except Exception, e:
        return None
register.filter('index_lookup', index_lookup)
register.filter('index', index_lookup)


@register.filter(name='in')
def in_(value, arg):
    """
    Filter that emulates "in" operator
    """
    return value in arg
register.filter('in', in_)


@register.filter(name='invalid')
def invalid(value, precision):
    """
    Filter that checks for invalid values 999.9 and 0 and replaces them with
    dashes.  valid values are formatted using floatformat and the inputed precision
    """
    if value:
        if value in (0,999.9):
            return '--'
        return floatformat(value, precision)
    return '--'
register.filter('invalid', invalid)


@register.filter(name='variable_dict_lookup')
def variable_dict_lookup(dict,key):
    """
    Allows using a variable for the dictionary and the key
    """
    if key in dict:
        return dict[key]
    else:
        return None
register.filter('variable_dict_lookup', variable_dict_lookup)


def format_atom(value):
    parts = []
    for i in value.split('_'):
        subscript = []
        for c in i[1:]:
            subscript.append(ROMAN_TO_GREEK[c] if c in ROMAN_TO_GREEK else c)
        parts.append('%s<sup>%s</sup>' % (i[0], ''.join(subscript)))
    return ''.join(parts)


@register.filter(name='sidechain_label')
def sidechain_label(value):
    """
    Filter that formats sidechain labels for angles and lengths
    """
    parts = [value[:3],': '] + [format_atom(value[5:])]
    return mark_safe(''.join(parts))
register.filter('sidechain_label', sidechain_label)


@register.filter(name='atom')
def atom(value):
    """
    Filter that formats atom bond lengths and angles
    """
    return mark_safe(format_atom(value))
register.filter('atom', atom)


@register.filter(name='label')
def label(value):
    """
    Formats any label converting from db field name to html label.
    sidechains are formated using sidechain_label().  if there is no
    matching label then the same value is returned
    """
    if value in PROPERTY_CHOICES_DICT:
        return PROPERTY_CHOICES_DICT[value]
    elif value[:10] == 'sidechain':
        return sidechain_label(value[10:])
    return value


@register.filter(name='sidechain_fields')
def sidechain_fields(form, key):
    """
    Filter that gets both the form field and hidden toggle field for the given key
    """
    fields = (str(form.__getitem__(key)), str(form.__getitem__(key+'_i')))
    return mark_safe(''.join(fields))
register.filter('sidechain_fields', sidechain_fields)


@register.filter(name='upper')
def upper(str):
    """
    Filter that converts string to uppercase
    """
    return str.upper()
register.filter('upper', upper)
