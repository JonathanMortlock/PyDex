"""PyDex types
Stefan Spence 13/04/20

Collection of functions for converting types
"""
import re
from distutils.util import strtobool

def BOOL(x):
    """Fix the conversion from string to Boolean.
    Any string with nonzero length evaluates as true 
    e.g. bool('False') is True. So we need strtobool."""
    try: return strtobool(x)
    except AttributeError: return bool(x)


def strlist(text):
    """Convert a string of a list of strings back into
    a list of strings."""
    return list(text[1:-1].replace("'","").split(', '))

def intstrlist(text):
    """Convert a string of a list of ints back into a list:
    (str) '[1, 2, 3]' -> (list) [1,2,3]"""
    try:
        return list(map(int, text[1:-1].split(',')))
    except ValueError: return []

def listlist(text):
    """Convert a string of nested lists into a
    list of lists."""
    return list(map(intstrlist, re.findall('\[[\d\s,]*\]', text)))

