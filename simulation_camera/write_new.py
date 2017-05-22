import sys
from PIL import Image

composite = str(sys.argv[1])

img = Image.new("RGBA",(2880,1980),(255, 255, 255, 255))
img.save(composite)
