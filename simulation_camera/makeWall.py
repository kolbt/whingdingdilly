import sys
from PIL import Image

composite = str(sys.argv[1])
dpi = int(sys.argv[2])

# Source images are 2000x2000
width = dpi * 16
height = dpi * 16

img = Image.new("RGBA", (width, height), (255, 255, 255, 255))
img.save(composite)

print("Creating image... ")
