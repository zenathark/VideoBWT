import numpy as np
import cv
import comp as cp
import wavelet as wv
import lwt

class codec(object):

    bpp = 2 ** 32
    dtype = int16
    decomp_level = 5

    def frame_compress(self, frame):
        comp_manager = cp.SPIHT(frame)
        comp_manager.bpp = self.bpp
        comp_manager.compress()
        c_frame = compressed_frame(comp_manager.output_stream, 
                comp_manager.n,
                (len(frame.data), len(frame.data)),
                frame.level)
        return c_frame 

    def frame_uncompress(self, c_frame, inplace = True):
        data = np.zeros(c_frame.size,dtype)
        wave = wv.wavelet2D(data,c_frame.decomp_level)
        comp_manager = cp.SPIHT(wave)
        if inplace:
            comp_manager.output_stream = c_frame.payload
        else:
            comp_manager.output_stream = list(c_frame)
        comp_manager.n = c_frame.max_n
        comp_manager.uncompress()
        return comp_manager

    def zip_from_cvimage(self, img):
        mtx = np.asarray(img)
        wave = lwt.cdf97(mtx)
        return frame_compress(wave)


    def zip_from_file(self, fname):
        img = cv.LoadImage(fname)
        return zip_from_image(img)

#TODO:
    def cframe_to_file(self):
        pass

    def unzip_from_file(self): #to a file
        pass

    def zip_from_video_file(self): #to a file
        pass

    def play_from_file(self)
        pass

class compressed_frame(object):
    def __init__(self, payload, max_n, size, decomp_level)
        self.payload = payload
        self.max_n = n
        self.size = size
        self.decomp_level = decomp_level

class video_capsule(object):
    buffer = []
    buffer_size = 5
    def __init__(self, max_n, size, decomp_level):

