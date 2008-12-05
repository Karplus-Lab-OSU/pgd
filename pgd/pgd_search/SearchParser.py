import re

"""
Validates a field to make sure that it has valid syntax for a search field
"""
def validateQueryField(str):
    r = re.compile(r'^(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?(,(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?)*$')
    m = r.match(str)
    return not m == None



