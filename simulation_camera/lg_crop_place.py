import sys
from PIL import Image

file = str(sys.argv[1])             # read in filename
comp = str(sys.argv[2])             # composite image name
loc_x = int(sys.argv[3])            # x-coord for paste
loc_y = int(sys.argv[4])            # y-coord for paste

img = Image.open(file)              # read in final tstep image
comp_img = Image.open(comp)

center = img.size[0] / 2
img2 = img.crop(
    (
        center - 900,
        center - 900,
        center + 900,
        center + 900,
    )
)                                   # crop away whitespace

#img3 = img2.resize((180,180))       # resize to prepare for paste
comp_img.paste(img2,(loc_x,loc_y))  # paste into composite image
comp_img.save(comp)                 # resave image
