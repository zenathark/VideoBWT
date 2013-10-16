import cv2
import os
import pickle


class VideoInfo:
    compression_algorithm = None
    frames = 0
    fps = 0


def split_raw(filename, dest_file):
    '''This function transcode a file to a Wavelet Compressed Video Format
    '''
    if dest_file[-1] != "/":
        dest_file += "/"
    dest_file += "raw/"
    original = cv2.VideoCapture(filename)
    loaded, frame = original.read()
    if not os.path.exists(dest_file):
        os.makedirs(dest_file)
    total_frames = original.get(cv2.cv.CV_CAP_PROP_FRAME_COUNT)
    current_frame = original.get(cv2.cv.CV_CAP_PROP_POS_FRAMES)
    while loaded and current_frame < total_frames:
        current_frame = original.get(cv2.cv.CV_CAP_PROP_POS_FRAMES)
        target_file = dest_file + str(int(current_frame)) + ".png"
        frame = cv2.cvtColor(frame, cv2.cv.CV_BGR2GRAY)
        cv2.imwrite(target_file, frame, [cv2.cv.CV_IMWRITE_PNG_COMPRESSION, 0])
        loaded, frame = original.read()
        current_frame = original.get(cv2.cv.CV_CAP_PROP_POS_FRAMES)
    info = VideoInfo()
    info.compression_algorithm = None
    info.frames = int(total_frames)
    info.fps = original.get(cv2.cv.CV_CAP_PROP_FPS)
    info_file = dest_file + "info.dat"
    pickle.dump(info, open(info_file, "wb"))


def compress_fvht(path, dest_path):
    print "TODO"


def compress_fvht_fullsearch(path, dest_path):
    print "TODO"


def compress_spiht(paath, dest_path):
    print "TODO"


def compress_spiht_fullsearch(path, dest_path):
    print "TODO"

if __name__ == '__main__':
    split_raw("/Users/juancgalan/video_test/akiyo_cif.mov",
              "/Users/juancgalan/Downloads/video_test/akiyo/")
