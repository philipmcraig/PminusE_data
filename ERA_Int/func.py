from __future__ import division
import pylab as pl
import os

def PrintFiles(path_to_files,start_of_file):
    """Function to isolate files that which start with a particular string
    within a particular directory then arranges them in the correct order.
    """
    list_of_files = os.listdir(path_to_files)
    filelist = []
    for each_file in list_of_files:
        if each_file.startswith(start_of_file):
            filelist.append(each_file)

    return filelist