import sys
from PIL import Image

composite = str(sys.argv[1])
factor = 5

# Make a 4x3 blank image
img = Image.new("RGBA",(2000 * factor, 1500 * factor),(255, 255, 255, 255))
img.save(composite)
