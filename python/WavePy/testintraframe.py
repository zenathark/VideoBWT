import intraframe as ifr
import numpy as np

frame = np.ones((20, 20)) * 2
keyframe = np.ones((20, 20)) * 2
frame[0:3, 0:3] = 1
keyframe[2:5, 2:5] = 1
error, m_vs = ifr.encode_motion_frame(frame, keyframe, (2, 2))
deframe = ifr.decode_motion_frame(error, m_vs, (2, 2), keyframe)
print deframe
