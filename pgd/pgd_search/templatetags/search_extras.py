from django.template.defaultfilters import stringfilter, floatformat
from django.utils.safestring import mark_safe
from django import template
register = template.Library()

from pgd_constants import AA_CHOICES_DICT


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


@register.filter(name='sidechain_label')
def sidechain_label(value):
    """
    Filter that formats sidechain labels for angles and lengths
    """
    parts = [value[:3],': '] + ['%s<sup>%s</sup>'%(i[0], i[1:]) for i in value[5:].split('_')]
    return mark_safe(''.join(parts))
register.filter('sidechain_label', sidechain_label)


@register.filter(name='atom')
def atom(value):
    """
    Filter that formats atom bond lengths and angles
    """
    parts = ['%s<sup>%s</sup>'%(i[0], i[1:]) for i in value.split('_')]
    return mark_safe(''.join(parts))
register.filter('atom', atom)


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
