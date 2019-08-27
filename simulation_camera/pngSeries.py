#import os
#from ovito.anim import *
#from ovito.data import *
#from ovito.io import import_file
#from ovito.modifiers import *
#from ovito.vis import *
#from ovito import dataset
#import ovito

# Read in a gsd, write out each frame as png

import os
import sys

import ovito
from ovito.io import *
from ovito.vis import *

infile = str(sys.argv[1])                               # gsd file
print(infile)
file_name = os.path.basename(infile)
outfile, file_extension = os.path.splitext(file_name)   # get base name

node = import_file(infile, multiple_frames = True)      # load file as trajectory
frames = node.source.num_frames
node.add_to_scene()

# Magenta
a_rcol = 1.0
a_gcol = 0.0
a_bcol = 1.0
a_alpha = 1.0

# Green
b_rcol = 0.0
b_gcol = 1.0
b_bcol = 0.0
b_alpha = 1.0

vp = ovito.dataset.viewports.active_vp


# Change this to .png
rs = RenderSettings(
    filename = outfile + "_frame_" + ".png",
    size = (2000,2000),
    background_color = (1.0, 1.0, 1.0),
    range = RenderSettings.Range.ANIMATION,
    renderer = OpenGLRenderer()
)                                                       # settings for render

vp.type = Viewport.Type.TOP                             # top view
vp.zoom_all()                                           # zoom to fit
vp.render(rs)                                           # render animation
