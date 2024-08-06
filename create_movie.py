##TODO:
##Options for more video formats

import os
import glob
import argparse

import cv2

DEFAULT_FPS = 10

VALID_FILE_EXTS = ['.png', '.PNG', '.jpg', '.jpeg', '.JPG', '.JPEG']

def create_movie_from_files(path_to_files, save_path, frames_per_second):
    """
    Create a movie from the desired list of images at each pressure level
    """
    # Each video has a frame per second which is number of frames in every second
    # Make list of files and sort
    files = []

    for ext in VALID_FILE_EXTS:
        print(f"Searching for files with extention {ext}")
        files = sorted(glob.glob(path_to_files + "*" + ext))
        if files:
            break

    if files:
        print("Files found.")
    else:
        raise ValueError(f"No files found for extensions {VALID_FILE_EXTS}")

    
    
    w, h = None, None
    for f in files:
        frame = cv2.imread(f)
    
        if w is None:
            # Setting up the video writer
            h, w, _ = frame.shape
            vid_spec = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
            writer = cv2.VideoWriter(save_path, vid_spec, frames_per_second, (w, h))
            
        writer.write(frame)
    
    writer.release()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=f'Convert sorted image files into a video. (Currently only MP4)')
    parser.add_argument('-t', '--title',
                        help='specifiy a title for the video file',
                        type=str,
                        default=None)
    parser.add_argument('--fps',
                        help='specifiy the frames per second for the video. Default = 10',
                        type=int,
                        default=None)
    parser.add_argument('input_file_directory',  
                        help='directory to read input files from',
                        type=str)
    parser.add_argument('save_directory',  
                        help='directory to save output video to',
                        type=str)

    args = parser.parse_args()

    if args.input_file_directory[-1:] != "/":
        input_dir = args.input_file_directory + "/"
    else:
        input_dir = args.input_file_directory

    print(input_dir)

    if args.save_directory[-1:] != "/":
        save_dir = args.save_directory + "/"
    else:
        save_dir = args.save_directory

    if args.title:
        save_path = os.path.join(save_dir, f'{args.title}.mp4')
    else:
        save_path = os.path.join(save_dir, f'Movie.mp4')

    print(save_path)


    if args.fps:
        fps = args.fps
    else:
        fps = DEFAULT_FPS

    print("Creating movie with specified settings...")
    create_movie_from_files(input_dir, save_path, fps)
    print("Done!")

