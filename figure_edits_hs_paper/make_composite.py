import sys
from PIL import Image

composite = 'SI_plane_figure.png'
nrows = 4
ncols = 8
width = int(sys.argv[1]) * nrows
height = int(sys.argv[2]) * ncols

# Make a blank image
img = Image.new("RGBA",(width, height),(255, 255, 255, 255))
img.save(composite)
