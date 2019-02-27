import sys
from PIL import Image

comp = str(sys.argv[1])
infile = str(sys.argv[2])
xpos = int(sys.argv[3])
ypos = int(sys.argv[4])

comp_img = Image.open(comp)         # open composite plane
img = Image.open(infile)            # open final simulation frame
comp_img.paste(img, (xpos, ypos))   # add final simulation frame to composite
comp_img.save(comp)                 # resave composite
