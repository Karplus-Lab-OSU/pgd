from django.template.defaultfilters import stringfilter
from django import template
register = template.Library()

from pgd_constants import AA_CHOICES_DICT

"""
Filter that returns the full AA code from the 1 letter AA code
"""
@register.filter(name='full_aa')
@stringfilter
def full_aa(value):
    return AA_CHOICES_DICT[value]
register.filter('full_aa',full_aa)