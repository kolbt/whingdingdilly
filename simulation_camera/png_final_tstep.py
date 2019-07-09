import os
import ovito
from ovito.anim import *
from ovito.data import *
from ovito.io import import_file
from ovito.modifiers import *
from ovito.vis import *

import sys

infile = str(sys.argv[1])
flip = int(sys.argv[2])
dpi = int(sys.argv[3])
outfile, file_extension = os.path.splitext(infile)      # get base name

node = import_file(infile, multiple_frames = True)      # load file as trajectory
node.add_to_scene()

a = 0
final = node.source.num_frames - 1                      # index of last frame

# Particle colors
r_tel = 0.294
g_tel = 0.654
b_tel = 0.611

r_orn = 0.807
g_orn = 0.650
b_orn = 0.321

# made this a while loop so that I could handle exceptions
while a == 0:
    # White background
    rs = RenderSettings(
        filename = outfile + ".png",
        size = (dpi, dpi),
        background_color = (1.0, 1.0, 1.0),
        renderer = OpenGLRenderer()
    )                                                       # settings for render

    ovito.dataset.anim.current_frame = final                # grab final snapshot
    vp = ovito.dataset.viewports.active_vp
    
    # Flip the color scheme A to orange, B to teal
    if flip:
        node.modifiers.append(SelectParticleTypeModifier(property='Particle Type',
                                                         types={0}))
        node.modifiers.append(AssignColorModifier(color=(r_tel, g_tel, b_tel)))

    # catch any corrupted data errors
    try:
        node.compute()                                          # set color of type A
        a = 1
    except Exception:
        final -= 1
        a = 0

    # Flip the color scheme A: teal, B: orange
    if flip:
        node.modifiers.append(SelectParticleTypeModifier(property='Particle Type',
                                                     types={1}))
        node.modifiers.append(AssignColorModifier(color=(r_orn, g_orn, b_orn)))

    # this block will only run if we haven't thrown an exception
    if a == 1:
        node.compute()                                          # set color of type B
        vp.type = Viewport.Type.TOP                             # top view
        vp.zoom_all()                                           # zoom to fit
        vp.render(rs)                                           # render image
