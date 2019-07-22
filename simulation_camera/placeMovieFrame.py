import sys
from PIL import Image

# Read in image and parameters
file = str(sys.argv[1])
comp = str(sys.argv[2])
pos = int(sys.argv[3])
factor = 4

# What are the coordinates? (2000x1500)
width = 2000 * factor
height = 1500 * factor
left = 0
midLeft = width / 4
midRight = width * 2 / 4
right = width * 3 / 4
top = 0
middle = height / 3
bottom = height * 2 / 3

# Assign x based on order in array
if pos == 0 or pos == 4 or pos == 8:
    x = left
elif pos == 1 or pos == 5 or pos == 9:
    x = midLeft
elif pos == 2 or pos == 6 or pos == 10:
    x = midRight
else:
    x = right
# Assign y based on order in array
if pos <= 3:
    y = top
elif 4 <= pos <= 7:
    y = middle
else:
    y = bottom


img = Image.open(file)              # read in final tstep image
comp_img = Image.open(comp)

img2 = img.resize((500*factor,500*factor))        # resize to prepare for paste
comp_img.paste(img2,(x,y))          # paste into composite image
comp_img.save(comp)                 # resave image
