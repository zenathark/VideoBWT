from scipy import *
from WavePy.tools import *
from WavePy.wavelet import *
import math
from numpy.numarray.numerictypes import Int

class speck(object):

    #wavelet object
    wv = 0

    #wavelet data
    dt = 0

    LIS = []
    LSP = []
    nextLIS = []
    S = []
    I = []
    n = 0
    output = []
    i_size_partition = 0
    bit_bucket = 0

    def __init__(self):
        pass
        
    def compress(self,wavelet, bpp):
        self.wv = wavelet
        self.dt = wavelet.data 
        self.bit_bucket = bpp * self.wv.rows * self.wv.cols
        self.initialization()
        wise_bit = self.n
        #sorting
        try:
            while self.n > 0:
                last_list = self.LIS.tail
                last_pixel = self.LSP.tail
                while self.LIS.index != last_list:
                    l = self.LIS.pop()
                    self.ProcessS(l)
                self.ProcessI()
                self.refinement(last_pixel)
                self.n -= 1
        except EOFError as e:
            print type(e)
            return [self.wv.cols, self.wv.rows, self.level, wise_bit, [self.output]]
        print "Elegant!!"
        return [self.wv.cols, self.wv.rows, self.level, wise_bit, [self.output]]
            
    def initialization(self):
        X = get_z_order(self.wv.rows * self.wv.cols)
        self.LIS = CircularStack(self.wv.cols * self.wv.rows)
        self.nextLIS = CircularStack(self.wv.cols * self.wv.rows)
        self.LSP = CircularStack(self.wv.cols * self.wv.rows)
        s_size = (self.wv.rows / 2**self.wv.level * self.wv.cols/ 2**self.wv.level)
        self.S = X[:s_size]
        del X[:s_size]
        self.I = X
        maxs = abs(self.wv.data)
        self.n = int(math.log(maxs.max(),2))
        self.LIS.push(self.S)
        self.i_partition_size = (self.wv.rows / 2**self.wv.level) ** 2

    def S_n(self,S):
        if len(S) == 0:
            return False
        T = np.array([i.tolist() for i in S])
        print len(S)
        return int((abs(self.dt[T[:,0],T[:,1]]).max() >= 2 ** self.n))

    def ProcessS(self,S):
        sn = self.S_n(S)
        self.out(sn)
        if sn == 1:
            if len(S) == 1:
                self.out(self.sign(S))
                self.push(S)
            else:
                self.CodeS(S)
        else:
            self.LIS.push(S)
    
    def CodeS(self,S):
        O = self.splitList(S)
        for o in O:
            sn = self.S_n(o)
            self.output += [sn]
            if sn == 1:
                if len(o) == 1:
                    self.out(self.sign(o))
                    self.push(o)
                else:
                    self.CodeS(o)
            else:
                self.LIS.push(o)
        pass

    def ProcessI(self):
        sn = self.S_n(self.I)
        if sn == 1:
            self.CodeI()

    def CodeI(self):
        part = self.splitList(self.I,self.i_partition_size)
        self.i_partition_size = self.i_partition_size * 4
        for i in range(3):
            self.ProcessS(part[i])
        self.I = part[3]
        self.ProcessI()
    
    def iInitialization(self, width, height, level, wise_bit):
        self.wv = wavelet2D(zeros((width,height), dtype=Int), level)
        self.dt = self.wv.data
        X = get_z_order(self.wv.rows * self.wv.cols)
        self.LIS = CircularStack(self.wv.cols * self.wv.rows)
        self.nextLIS = CircularStack(self.wv.cols * self.wv.rows)
        self.LSP = CircularStack(self.wv.cols * self.wv.rows)
        s_size = (self.wv.rows / 2**self.wv.level * self.wv.cols/ 2**self.wv.level)
        self.S = X[:s_size]
        del X[:s_size]
        self.I = X
        self.n = wise_bit
        self.LIS.push(self.S)
        self.i_partition_size = (self.wv.rows / 2**self.wv.level) ** 2
    
    def expand(self,stream, width, height, level, wise_bit):
        self.iInitialization(width, height, level, wise_bit)
        self.output = stream
        #sorting
        try:
            while self.n > 0:
                last_list = self.LIS.tail
                last_pixel = self.LSP.tail
                while self.LIS.index != last_list:
                    l = self.LIS.pop()
                    self.iProcessS(l)
                self.iProcessI()
                self.iRefinement(last_pixel)
                self.n -= 1
        except EOFError as e:
            print type(e)
            return self.wv
        print "Elegant!!"
        return self.wv
    
    def iProcessS(self,S):
        sn = self.read()
        if sn == 1:
            if len(S) == 1:
                self.dt[S[0][0],S[0][1]] = 2 ** self.n * -1 * self.read() 
                self.push(S)
            else:
                self.iCodeS(S)
        else:
            self.LIS.push(S)
    
    def iCodeS(self,S):
        O = self.splitList(S)
        for o in O:
            sn = self.read()
            if sn == 1:
                if len(o) == 1:
                    self.dt[S[0][0],S[0][1]] = 2 ** self.n * -1 * self.read() 
                    self.push(o)
                else:
                    self.iCodeS(o)
            else:
                self.LIS.push(o)
        pass
            
    def iProcessI(self):
        sn = self.read()
        if sn == 1:
            self.iCodeI()
            
    def iCodeI(self):
        part = self.splitList(self.I,self.i_partition_size)
        self.i_partition_size = self.i_partition_size * 4
        for i in range(3):
            self.iProcessS(part[i])
        self.I = part[3]
        self.iProcessI()
    
    def sign(self,S):
        if S[0] >= 0:
            return 0
        else:
            return 1
    
    def splitList(self,l, size = 0):
        if size == 0:
            if len(l) % 4 != 0:
                raise IndexError 
            size = int(len(l) / 4)
            return [l[i*size:(i+1)*size] for i in (0,1,2,3)]
        else:
            return [l[i*size:(i+1)*size] for i in (0,1,2)] + [l[:-size*3]]
        
    def out(self,data):
        if len(self.output) < self.bit_bucket:
            self.output += [data]
        else:
            raise EOFError
    
    def read(self):
        if len(self.output) > 0:
            return self.output.pop(0)
        else:
            raise EOFError
            
    def refinement(self, end):
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            if (self.dt[i[0],i[1]] & self.n) > 0:
                self.out(1)
            else:
                self.out(0)
            c = (c + 1) % self.LSP.size
            
    def iRefinement(self, end):
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            if (self.read()) > 0:
                self.dt[i[0],i[1]] |= 2 ** self.n
            c = (c + 1) % self.LSP.size
                
    def push(self, data):
        self.LSP.push(data[0])
        
class fv_speck(speck):
    
    def __init__(self):
        pass
    
    def compress(self,wavelet, bpp, fp, window):
        self.calculate_fovea_length()
        super(fv_speck,self).compress(self.wv, bpp)
        
    
    def push(self, data):
        fv = self.calculate_fovea_w(data)
        if fv >= self.get_current_bpp():
            self.LSP.push(data)
            
    def refinement(self):
        for i in self.LSP:
            fv = self.calculate_fovea_w(i)
            if fv >= self.get_current_bpp():
                if (self.dt[i[0],i[1]] & self.n) > 0:
                    self.out(1)
                else:
                    self.out(0)
        
    def calculate_fovea_w(self, ij):
        try:
            P = self.get_center(ij)
        except NameError:
            return self.Lbpp
        d = self.norm(P[1] - ij[1],P[0] - ij[0]) * 2**P[2] / self.fovea_length
        if d<self.alpha:
            return self.Lbpp
        elif d>=1:
            return self.lbpp
        else:
            return self.powerlaw(d) * (self.Lbpp - self.lbpp) + self.lbpp
    
    def get_center(self,ij):
        if (ij[0] == 0 and ij[1] == 0):
            raise NameError("ij on HH")
        else:
            if ij[0] == 0:
                aprx_level_r = self.wavelet.level + 1
            else:
                aprx_level_r = math.ceil(math.log(self.wavelet.rows/float(ij[0]),2))
                if aprx_level_r > self.wavelet.level:
                    aprx_level_r = self.wavelet.level + 1
            if ij[1] == 0:
                aprx_level_c = self.wavelet.level + 1
            else:
                aprx_level_c = math.ceil(math.log(self.wavelet.rows/float(ij[1]),2))
                if aprx_level_c > self.wavelet.level:
                    aprx_level_c = self.wavelet.level + 1
            if (aprx_level_r > self.wavelet.level) and (aprx_level_c > self.wavelet.level):
#                raise NameError("ij on HH")
                y = float(self.P[0]) / 2**(aprx_level_r-1)
                x = float(self.P[1]) / 2**(aprx_level_r-1)
                return (y,x,aprx_level_r-1) 
            if aprx_level_r <= aprx_level_c:
                aprx_level = aprx_level_r
            else:
                aprx_level = aprx_level_c
            y = float(self.P[0]) / 2**aprx_level
            x = float(self.P[1]) / 2**aprx_level 
            if aprx_level_r == aprx_level:
                y += float(self.wavelet.rows) / 2**aprx_level
            if aprx_level_c == aprx_level:
                x += float(self.wavelet.cols) / 2**aprx_level
            return (y, x, aprx_level)
        
    def calculate_fovea_length(self):
        H = len(self.wavelet.data)
        W = len(self.wavelet.data[0])
        k = zeros(4)
        k[0] = self.norm(self.P[0],H-self.P[1])
        k[1] = self.norm(W-self.P[0],self.P[1])
        k[2] = self.norm(W-self.P[0],H-self.P[1])
        k[3] = self.norm(self.P[0],H-self.P[1])
        self.fovea_length = k.max()
        
    def printFoveaWindow(self):
        window = zeros((self.wavelet.rows, self.wavelet.cols))
        points = self.wv.get_z_order(self.wavelet.rows * self.wavelet.cols)
        for i in points:
            window[tuple(i)] = self.calculate_fovea_w(i)
        return window
    
    def get_current_bpp(self):
        bpp = len(self.output)
        bpp /= float(self.wavelet.rows*self.wavelet.cols)
        return bpp