#import os
#from ovito.anim import *
#from ovito.data import *
#from ovito.io import import_file
#from ovito.modifiers import *
#from ovito.vis import *
#from ovito import dataset
#import ovito

import os
import sys

import ovito
from ovito.io import *
from ovito.vis import *

infile = str(sys.argv[1])                               # gsd file
file_name = os.path.basename(infile)
outfile, file_extension = os.path.splitext(file_name)   # get base name

print(outfile)

node = import_file(infile, multiple_frames = True)      # load file as trajectory
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

rs = RenderSettings(
    filename = outfile + ".mp4",
    size = (1000,1000),
    background_color = (1.0, 1.0, 1.0),
    range = RenderSettings.Range.ANIMATION,
    renderer = OpenGLRenderer()
)                                                       # settings for render

#node.modifiers.append(SelectParticleTypeModifier(property='Particle Type',
#                                                 types={0}))
#node.modifiers.append(AssignColorModifier(color=(a_rcol, a_gcol, a_bcol)))
#
#node.compute()                                          # set color of type A
#
#node.modifiers.append(SelectParticleTypeModifier(property='Particle Type',
#                                                 types={1}))
#node.modifiers.append(AssignColorModifier(color=(b_rcol, b_gcol, b_bcol)))
#node.compute()                                          # set color of type B

vp.type = Viewport.Type.TOP                             # top view
vp.zoom_all()                                           # zoom to fit
vp.render(rs)                                           # render animation
