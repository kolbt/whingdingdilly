import sys
from PIL import Image

composite = str(sys.argv[1])

img = Image.new("RGBA",(28800,19800),(255, 255, 255, 255))
#img = Image.new("RGBA",(37800,1800),(255, 255, 255, 255))
img.save(composite)
