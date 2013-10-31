import tools as ts
import wavelet as wvt
import math
from numpy.numarray.numerictypes import Int
import logging
import numpy as np


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
    log = []
    debug = True
    f_log = []
    logger = None
    _idx = 0
    _logidx = 0
    out_idx = []

    def __init__(self):
        pass

    def compress(self, wavelet, bpp):
        self.wv = wavelet
        self.dt = wavelet.data
        self.bit_bucket = bpp * self.wv.rows * self.wv.cols
        self.initialization()
        wise_bit = self.n
        #sorting
        try:
            while self.n > 0:
                print self.n
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
            return [self.wv.cols, self.wv.rows, self.wv.level,
                    wise_bit, self.output]
        #print "Elegant!!"
        return [self.wv.cols, self.wv.rows, self.wv.level,
                wise_bit, self.output]

    def initialization(self):
        X = wvt.get_z_order(self.wv.rows * self.wv.cols)
        self.LIS = ts.CircularStack(self.wv.cols * self.wv.rows)
        self.nextLIS = ts.CircularStack(self.wv.cols * self.wv.rows)
        self.LSP = ts.CircularStack(self.wv.cols * self.wv.rows)
        s_size = (self.wv.rows * self.wv.cols / 2 ** (2 * self.wv.level))
        self.S = X[:s_size]
        del X[:s_size]
        self.I = X
        maxs = abs(self.wv.data)
        self.n = int(math.log(maxs.max(), 2))
        self.LIS.push(self.S)
        self.i_partition_size = (self.wv.rows / 2 ** self.wv.level) ** 2
        self.output = [0] * self.bit_bucket
        self.out_idx = 0

    def S_n(self, S):
        if len(S) == 0:
            return False
        T = np.array([i.tolist() for i in S])
        return int((abs(self.dt[T[:, 0], T[:, 1]]).max() >= 2 ** self.n))

    def ProcessS(self, S):
        sn = self.S_n(S)
        self.out(sn)
        #self.writeLog("ProcessS", "Sn", "S", len(S), sn)
        if sn == 1:
            if len(S) == 1:
                self.out(self.sign(S))
                self.push(S)
            else:
                self.CodeS(S)
        else:
            self.LIS.push(S)

    def CodeS(self, S):
        O = self.splitList(S)
        for o in O:
            sn = self.S_n(o)
            self.out(sn)
            #self.writeLog("CodeS", "Sn", "S", len(S), sn)
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
        self.out(sn)
        if sn == 1:
            self.CodeI()

    def CodeI(self):
        part = self.splitList(self.I, self.i_partition_size)
        self.i_partition_size = self.i_partition_size * 4
        for i in range(3):
            self.ProcessS(part[i])
        self.I = part[3]
        self.ProcessI()

    def iInitialization(self, width, height, level, wise_bit):
        self.wv = wvt.wavelet2D(np.zeros((width, height), dtype=Int), level)
        self.dt = self.wv.data
        self.wv.level = level
        X = wvt.get_z_order(self.wv.rows * self.wv.cols)
        self.LIS = ts.CircularStack(self.wv.cols * self.wv.rows)
        self.nextLIS = ts.CircularStack(self.wv.cols * self.wv.rows)
        self.LSP = ts.CircularStack(self.wv.cols * self.wv.rows)
        s_size = (self.wv.rows * self.wv.cols / 2 ** (2 * self.wv.level))
        self.S = X[:s_size]
        del X[:s_size]
        self.I = X
        self.n = wise_bit
        self.LIS.push(self.S)
        self.i_partition_size = (self.wv.rows / 2 ** self.wv.level) ** 2
        self._idx = 0
        self._logidx = 0

    def expand(self, stream, width, height, level, wise_bit):
        self.iInitialization(width, height, level, wise_bit)
        self.output = stream
        #sorting
        try:
            while self.n > 0:
                print self.n
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
        #print "Elegant!!"
        return self.wv

    def iProcessS(self, S):
        sn = self.read()
        if sn == 1:
            if len(S) == 1:
                sg = self.read()
                self.createCoeff(S[0], sg)
                self.push(S)
            else:
                self.iCodeS(S)
        else:
            self.LIS.push(S)

    def iCodeS(self, S):
        O = self.splitList(S)
        for o in O:
            sn = self.read()
            if sn == 1:
                if len(o) == 1:
                    sg = self.read()
                    self.createCoeff(o[0], sg)
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
        part = self.splitList(self.I, self.i_partition_size)
        self.i_partition_size = self.i_partition_size * 4
        for i in range(3):
            self.iProcessS(part[i])
        self.I = part[3]
        self.iProcessI()

    def sign(self, S):
        if self.dt[S[0][0], S[0][1]] >= 0:
            return 0
        else:
            return 1

    def splitList(self, l, size=0):
        if size == 0:
            if len(l) % 4 != 0:
                raise IndexError
            size = int(len(l) / 4)
            return [l[i * size:(i + 1) * size] for i in (0, 1, 2, 3)]
        else:
            return [l[i * size:(i + 1) * size] for i in (0, 1, 2)] \
                + [l[size * 3:]]

    def out(self, data):
        if self.out_idx < self.bit_bucket:
            self.output[self.out_idx] = data
            self.out_idx += 1
        else:
            raise EOFError

    def read(self):
        if self._idx < len(self.output):
            self._idx += 1
            return self.output[self._idx - 1]
        else:
            raise EOFError

    def refinement(self, end):
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            if self.dt[i[0], i[1]] > 0:
                coeff = self.dt[i[0], i[1]]
            else:
                coeff = abs(self.dt[i[0], i[1]])
            if (coeff & 2 ** self.n) > 0:
                self.out(1)
            else:
                self.out(0)
            c = (c + 1) % self.LSP.size

    def iRefinement(self, end):
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            if (self.read()) > 0:
                if self.dt[i[0], i[1]] > 0:
                    self.dt[i[0], i[1]] |= 2 ** self.n
                else:
                    self.dt[i[0], i[1]] = (abs(self.dt[i[0], i[1]])
                                           | 2 ** self.n) * -1
            c = (c + 1) % self.LSP.size

    def push(self, data):
        self.LSP.push(data[0])

    def createCoeff(self, coords, sg, wv=None):
        if wv is None:
            self.dt[coords[0], coords[1]] = (2 ** self.n) * \
                ((sg * 2 - 1) * -1)

    def writeLog(self, method, reason, obj, size, value):
        if self.debug:
            self.log += [method, reason, obj, size, value]


class fv_speck(speck):

    def __init__(self):
        pass

    def compress(self, wavelet, bpp, lbpp, f_center, alpha, c, gamma):
        self.Lbpp = bpp
        self.lbpp = lbpp
        self.alpha = alpha
        self.wv = wavelet
        self.dt = wavelet.data
        self.P = f_center
        self.c = c
        self.gamma = gamma
        self.calculate_fovea_length()
        return super(fv_speck, self).compress(self.wv, bpp)

    def expand(self, stream, width, height, level, wise_bit, bpp, lbpp,
               f_center, alpha, c, gamma):
        self.Lbpp = bpp
        self.lbpp = lbpp
        self.alpha = alpha
        self.P = f_center
        self.c = c
        self.wv = wvt.wavelet2D(np.zeros((width, height), dtype=Int), level)
        self.dt = self.wv.data
        self.wv.level = level
        self.gamma = gamma
        self.calculate_fovea_length()
        return super(fv_speck, self).expand(stream, width, height, level,
                                            wise_bit)

    def refinement(self, end):
        print('iRefinement I' + str(len(self.I)) + ' ' + str(len(self.output)))
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            fv = self.calculate_fovea_w(i)
            if fv >= self.get_current_bpp():
                if self.dt[i[0], i[1]] > 0:
                    coeff = self.dt[i[0], i[1]]
                else:
                    coeff = abs(self.dt[i[0], i[1]])
                if (coeff & 2 ** self.n) > 0:
                    self.out(1)
                else:
                    self.out(0)
            c = (c + 1) % self.LSP.size

    def iRefinement(self, end):
        c = self.LSP.index
        while c != end:
            i = self.LSP.data[c]
            fv = self.calculate_fovea_w(i)
            if fv >= (self.get_dec_bpp()):
                if (self.read()) > 0:
                    if self.dt[i[0], i[1]] > 0:
                        self.dt[i[0], i[1]] |= 2 ** self.n
                    else:
                        self.dt[i[0], i[1]] = (abs(self.dt[i[0], i[1]]) |
                                               2 ** self.n) * -1
            c = (c + 1) % self.LSP.size

    def calculate_fovea_w(self, ij):
        try:
            P = self.get_center(ij)
        except NameError:
            return self.Lbpp
        d = self.norm(P[1] - ij[1], P[0] - ij[0]) * 2 ** P[2] / \
            self.fovea_length
        if d < self.alpha:
            return self.Lbpp
        elif d >= 1:
            return self.lbpp
        else:
            return self.powerlaw(d) * (self.Lbpp - self.lbpp) + self.lbpp

    def get_center(self, ij):
        if (ij[0] == 0 and ij[1] == 0):
            raise NameError("ij on HH")
        else:
            if ij[0] == 0:
                aprx_level_r = self.wv.level + 1
            else:
                aprx_level_r = math.ceil(math.log(self.wv.rows /
                                                  float(ij[0]), 2))
                if aprx_level_r > self.wv.level:
                    aprx_level_r = self.wv.level + 1
            if ij[1] == 0:
                aprx_level_c = self.wv.level + 1
            else:
                aprx_level_c = math.ceil(math.log(self.wv.rows /
                                                  float(ij[1]), 2))
                if aprx_level_c > self.wv.level:
                    aprx_level_c = self.wv.level + 1
            if (aprx_level_r > self.wv.level) and \
               (aprx_level_c > self.wv.level):
#                raise NameError("ij on HH")
                y = float(self.P[0]) / 2 ** (aprx_level_r - 1)
                x = float(self.P[1]) / 2 ** (aprx_level_r - 1)
                return (y, x, aprx_level_r - 1)
            if aprx_level_r <= aprx_level_c:
                aprx_level = aprx_level_r
            else:
                aprx_level = aprx_level_c
            y = float(self.P[0]) / 2 ** aprx_level
            x = float(self.P[1]) / 2 ** aprx_level
            if aprx_level_r == aprx_level:
                y += float(self.wv.rows) / 2 ** aprx_level
            if aprx_level_c == aprx_level:
                x += float(self.wv.cols) / 2 ** aprx_level
            return (y, x, aprx_level)

    def calculate_fovea_length(self):
        H = self.wv.rows
        W = self.wv.cols
        k = np.zeros(4)
        k[0] = self.norm(self.P[0], H - self.P[1])
        k[1] = self.norm(W - self.P[0], self.P[1])
        k[2] = self.norm(W - self.P[0], H - self.P[1])
        k[3] = self.norm(self.P[0], H - self.P[1])
        self.fovea_length = k.max()

    def printFoveaWindow(self):
        window = np.zeros((self.wv.rows, self.wv.cols))
        points = self.wv.get_z_order(self.wavelet.rows * self.wavelet.cols)
        for i in points:
            window[tuple(i)] = self.calculate_fovea_w(i)
        return window

    def get_current_bpp(self):
        bpp = len(self.output)
        bpp /= float(self.wv.rows * self.wv.cols)
        return bpp

    def get_dec_bpp(self):
        bpp = self._idx
        bpp /= float(self.wv.rows * self.wv.cols)
        return bpp

    def norm(self, x, y):
        mx = abs(x)
        if mx < abs(y):
            mx = abs(y)
        return mx  # math.sqrt(float(x**2 + y ** 2))

    def powerlaw(self, n):
        return self.c * (1 - ((n - self.alpha) / (1 - self.alpha))) \
            ** self.gamma
