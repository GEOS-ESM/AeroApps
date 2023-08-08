import os
try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk
import fnmatch

def find(path=os.getcwd(), ext='.txt'):
    '''Recursive search function top-down.'''
    for (root, dirs, files) in walk(path):
        for f in fnmatch.filter(files, '*'+ext):
            yield os.path.join(root, f)

def find_fast(path=os.getcwd(), ext='.txt'):
    '''Recursive search function top-down.'''
    for (root, dirs, files) in walk(path):
        for d in dirs:
            yield os.path.join(root, d)
        #for f in fnmatch.filter(files, '*'+ext):
        #    yield os.path.join(root, f)
