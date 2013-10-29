"""
..  module::intraframe
    :platform: Unix
    :synopsis: This module has metods for intraframe encoding of video

..  moduleauthor:: Juan Carlos Galan Hernandez <juan.galanhz@udlap.mx>

"""
import itertools as it
import numpy as np


def encode_motion_frame(frame, key_frame, window_size, region='FULL'):
    """this function encode a video frame using motion vectors.

    Args:
        frame (np.array): A two dimentional matrix that represents the frame
        to be encoded.
        key_frame: The reference frame used to encode the given frame
        window_size (tuple): The size of the size of the macroblock to be
        used for encoding given as (row size, column size)
        region (tuple): The size of the searchin region given as (row size,

    Kwargs:
        column size). Instead, if 'FULL' is given, the size of the searching
        region will be the full key_frame

    Returns:
        The motion error frame matrix and an array of motion vector returned as
        a tuple of the form:
            (np.array, [(tuple, tuple, ...)])
    """
    #frame size
    k_r, k_c = key_frame.shape
    #window size
    r = c = window_size
    encoded_frame = np.zeros(frame.shape)
    motion_vectors = []
    for i, j in it.product(range(0, k_r, r), range(0, k_c, c)):
        N = i + r
        M = j + c
        encoded_frame[i:N, j:M], motion = \
            fullsearch(frame[i:N, j:M], (i, j), key_frame, region)
        motion_vectors += [motion]
    return (encoded_frame, motion_vectors)


def fullsearch(window, window_coord, key_frame, region="FULL"):
    """This function calculates the motion estimation vector of the given
    vector over the key frame using full search.

    Args:
        window (np.array): A two dimentional matrix used for calculating the
        motion vector
        win_coord (tuple): A tuple representing the position of the window in
        its original frame
        key_frame (np.array): A two dimentional matrix that holds the key frame

    Kwargs:
        region: The column and row length of the region to be searched.
        The string 'FULL' causes the function to search over all the key frame.
        If an int is given, the function will delimite the search over an
        squared area with radius of region centered at win_coords.

    Returns:
        Returns a tuple (e, (u,v)) that is the best error prediction 'e' found
        and its motion vector where the given window with coordinates (u,v).
    """
    y, x = window_coord
    #row size of the window
    r, c = window.shape
    if region == "FULL":
        y_0 = 0
        x_0 = 0
        y_t = key_frame.shape[0]
        x_t = key_frame.shape[1]
    else:
        y_r = x_r = region
        #check if region search is not outside of the key_frame
        y_0 = y - y_r if (y - y_r) >= 0 else 0
        x_0 = x - x_r if (x - x_r) >= 0 else 0
        y_t = y + y_r if (y + y_r) <= key_frame.shape[0] \
            else key_frame.shape[0]
        x_t = x + x_r if (x + x_r) <= key_frame.shape[1] \
            else key_frame.shape[1]
    #error stores the motion estimation error calculated as a difference of the
    #currento window minus an area of the key_frame
    mae, error = MAE(window, key_frame[y:y + r, x:x + c])
    best_prediction_coord = (y, x)
    for i, j in it.product(range(y_0, y_t, window.shape[0]),
                           range(x_0, x_t, window.shape[1])):
        local_mae, local_error = MAE(window, key_frame[i:i + r, j:j + c])
        if local_mae < mae:
            mae = local_mae
            best_prediction_coord = (i, j)
            error = local_error
    best_y = best_prediction_coord[0] - y
    best_x = best_prediction_coord[1] - x
    return (error, (best_y, best_x))


def MAE(w, kw):
    """This function calculates the Mean Absolute Error over two given matrix.
    The used formula is defined as

        MAE(w) = 1/(n^2) \sum{k=0}{N-1}\sum{i=0}{N-1}\abs{w[k,i]-kw[k,i]}

    where w is the first matrix (the window used for calculating the motion
    vector) and kw is the extracted window from the key frame.

    Args:
        w (np.array): A two dimensional matrix
        kw (np.array): Another two dimensional matrix

    Returns:
        This function returns a tuple that contains a float as a result of
        the MAE algorithm as a first element and a 2D matrix of the motion
        estimation error
    """
    N = w.shape[0]
    error = w - kw
    mae = 1.0 / N ** 2 * sum(sum(abs(w - kw))) * 2
    return (mae, error)


def decode_motion_frame(error, motion_vectors, window_size, key_frame):
    '''This method decode an encoded frame using motion estimation
    '''
    #frame size
    k_r, k_c = key_frame.shape
    #window size
    r = c = window_size
    decoded_frame = np.zeros(error.shape)
    #composed index for easy iteration
    index = zip(it.product(range(0, k_r, r), range(0, k_c, c)), motion_vectors)
    for (i, j), (d_i, d_j) in index:
        N = i + r
        M = j + r
        window = error[i:N, j:M]
        key_window = key_frame[i + d_i:N + d_i, j + d_j:M + d_j]
        decoded_frame[i:N, j:M] = window + key_window
    return decoded_frame
