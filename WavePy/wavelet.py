'''
Created on 02/06/2012

@author: zenathar
'''
from numpy import array 

class wavelet(object):
    def __init__(self):
        print("echo")

class wt2D(object):
    '''
    classdocs
    '''
    
    level = 1
    w_type = float


    def __init__(self, data, level, w_type):
        '''
        Constructor
        '''
        if not isinstance(data, array) or len(data) != 2:
            raise TypeError, "data must be a 2D numpy array matrix"
        self.level = level
        self.type = w_type
        self.data = data
        