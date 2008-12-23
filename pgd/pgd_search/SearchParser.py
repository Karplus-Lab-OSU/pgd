import re

"""
Validates a field to make sure that it has valid syntax for a search field
"""
def validateQueryField(str):
    return re.compile(r'^(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?(,(-?([1-9]\d*|0)(\.\d+)?|(\.\d+))(-(-?([1-9]\d*|0)(\.\d+)?|(\.\d+)))?)*$').match(str) != None
