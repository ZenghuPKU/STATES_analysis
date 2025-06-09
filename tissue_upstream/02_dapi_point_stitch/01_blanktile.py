import numpy as np
from tifffile import imwrite
import os

user_dir = "/home/huzeng_pkuhpc/gpfs3/yly/20250114mousebrain/alldata"
sample = "B4_WT_mousebrain2"
subdirs = {
    "source_data": "01_data",
    "registration": "02_registration",
    "segmentation": "03_segmentation",
    "stitch": "04_stitch",
    "final": "05_results"
}

def getOutpath(subdir):
    path = os.path.join(user_dir, sample, subdirs[subdir])
    return path
    
blank = np.zeros((3072, 3072), dtype=np.int8)

# Create output directory if it doesn't exist
outpath = getOutpath('stitch')
os.makedirs(outpath, exist_ok=True)

# Save blank tile to output directory
imwrite(f"{outpath}/blank.tif", blank)
