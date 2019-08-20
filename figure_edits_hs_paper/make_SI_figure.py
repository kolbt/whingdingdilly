import sys
from PIL import Image

# Paths of raw vs scatter
composite = '/Users/kolbt/Desktop/ALL_FIGURES/supplemental_info/all_phase_planes/SI_plane_figure.png'
rawPth = '/Users/kolbt/Desktop/ALL_FIGURES/supplemental_info/all_phase_planes/sim_overlay/'
sctPth = '/Users/kolbt/Desktop/ALL_FIGURES/supplemental_info/all_phase_planes/mips_algorithm/'

# Get PeB plane
pb = int(sys.argv[1])
file = 'HS_peB_' + str(pb) + '.png'
# All images are 6400 x 4800
width = int(sys.argv[2])
height = int(sys.argv[3])

# Read in the raw file
img1 = Image.open(rawPth + file)
img1 = img1.resize((width,height))
# Read in the scatter file
img2 = Image.open(sctPth + file)
img2 = img2.resize((width,height))
# Open the composite
comp = Image.open(composite)

if pb % 20 == 0:
    x1 = 0
else:
    x1 = 2 * width

y1 = int(pb / 20) * height
x2 = x1 + width
y2 = y1

# Crop and paste to composite file
comp.paste(img1,(x1,y1))
comp.paste(img2,(x2,y2))
comp.save(composite)
