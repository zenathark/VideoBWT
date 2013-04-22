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

    def __init__(self, wavelet):
        wv = wavelet
        dt = wavelet.data
        LIS = CircularStack(wv.cols * wv.rows)
        LSP = CircularStack(wv.cols * wv.rows)
        
    def initialization(self):
        X = self.wv.get_morton_order(self.wv.rows * self.wv.cols)
        self.LSP = []
        s_size = (self.wv.row / 2**wv.level * self.wv.cols/ 2**wv.level)
        self.S = X[:s_size-1]
        del X[:s_size-1]
        I = X
        maxs = abs(self.wavelet.data)
        self.n = int(math.log(maxs.max(),2))
        self.LIS = []
        self.LIS += [(self.S)]


    def compress(self):


