import sys
from PIL import Image

composite = str(sys.argv[1])

# Source images are 2000x2000
width = 250 * 16
height = 250 * 11

img = Image.new("RGBA", (width, height), (255, 255, 255, 255))
img.save(composite)

print("Creating image... ")
