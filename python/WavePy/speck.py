from scipy import *
from WavePy.tools import *

class speck(object):

    #wavelet object
    wv = 0

    #wavelet data
    dt = 0

    LIS = []
    LSP = []
    S = []
    I = []
    n = 0
    output = []

    def __init__(self, wavelet):
        self.wv = wavelet
        self.dt = wavelet.data
        
    def initialization(self):
        X = self.wv.get_morton_order(self.wv.rows * self.wv.cols)
        self.LIS = CircularStack(self.wv.cols * self.wv.rows)
        self.LSP = CircularStack(self.wv.cols * self.wv.rows)
        s_size = (self.wv.row / 2**self.wv.level * self.wv.cols/ 2**self.wv.level)
        self.S = X[:s_size-1]
        del X[:s_size-1]
        self.I = X
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        self.LIS.push(self.S)

    def S_n(self,S):
        pass

    def ProcessS(self,S):
        sn = self.S_n(S)
        self.output += [sn]
        if sn == 1:
            if len(S) == 1:
                self.output += [sign(S)]
                self.LSP.push(S)
                #TODO erase S
            else:
                self.CodeS(S)
        #TODO add S
    
    def CodeS(self,S):
        pass

    def ProcessI(self):
        pass

    def CodeI(self):
        pass
    
    def sign(self,S):
        pass