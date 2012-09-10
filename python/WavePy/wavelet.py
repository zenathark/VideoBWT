'''
Created on 02/06/2012

@author: zenathar
'''
import numpy as np
from numpy import array 
import pdb

class wavelet(object):
    def __init__(self):
        print("echo")

class wavelet2D(object):
    '''
    classdocs
    '''
    
    level = 1
    w_type = float


    def __init__(self, data, level, w_type = float):
        '''
        Constructor
        '''
        if not isinstance(data, np.ndarray) or data.ndim != 2:
            raise TypeError, "data must be a 2D numpy array matrix"
        self.level = level
        self.w_type = w_type
        self.data = data
        
    def getHH(self):
        '''
        Returns the HH band of the wavelet
        '''
        data1 = self.data[0:len(self.data) / 2**self.level,0:len(self.data[0]) / 2 ** self.level]
        return data1

    def getLH(self,subdata):
        '''
        Returns the subband LH of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[len(subdata) / 2:,:len(subdata[0]) / 2]
        return data1

    def getHL(self,subdata):
        '''
        Returns the subband HL of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[:len(subdata) / 2,len(subdata[0]) / 2:]
        return data1

    def getLL(self,subdata):
        '''
        Returns the subband LH of a wavelet assuming one level of decomposition
        '''
        data1 = subdata[len(subdata) / 2:,len(subdata[0]) / 2:]
        return data1

    def getLH_n(self, level):
        '''
        Returns the LH subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getLH(self.data[:rows, :cols])
        return _subband
    
    def getHL_n(self, level):
        '''
        Returns the HL subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getHL(self.data[:rows, :cols])
        return _subband

    def getLL_n(self, level):
        '''
        Returns the LL subband of level "level" of the wavelet decomposition
        '''
        level = level - 1
        #calculate the number of rows and columns that the decomposition must have at level - 1
        rows = len(self.data) / 2 ** level
        cols = len(self.data[0]) / 2 ** level
        _subband = self.getLL(self.data[:rows, :cols])
        return _subband

