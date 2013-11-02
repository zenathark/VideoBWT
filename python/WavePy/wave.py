""" Wavelet definition.

File: wave.py
Author: jcgalanh@gmail.com
Description: Definition of the wavelet data type

"""

import numpy as np
import WavePy.tools as tools


class Wavelet(np.ndarray):

    """
    This object represents a wavelet.

    The fundamental element for signal processing using wavelets is an N matrix
    that holds the coefficients of a wavelet decomposition. This object extends
    from numpy.ndarray and extends it to hold the extra values needed for a
    wavelet data set

    """

    def __init__(self, level, **kwargs):
        if 'filter' in kwargs:
            filter_ = kwargs['filter']
        else:
            filter_ = None
        np.ndarray.__init__(self, kwargs)


    @staticmethod
    def fromArray(array, level, filter_=None):
        """Create a wavelet.

        This method creates a wavelet object using a numpy.ndarray as base

        Args:
            array. A numpy.ndarray as a base for this wavelet
            level. Level of decomposition of this wavelet
            filter. Filter bank name used

        Return:
            A Wavelet object with the same data as the numpy.ndarray object.
            The data is shared between both objects

        """
        tools.check_ndarray(array)
        wavelet = Wavelet(level, filter=filter_, shape=array.shape,
                          buffer=array, dtype=array.dtype)
        return wavelet
