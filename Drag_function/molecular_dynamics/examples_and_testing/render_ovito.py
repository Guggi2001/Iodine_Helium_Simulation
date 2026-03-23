from ovito.io import import_file
from ovito.modifiers import AssignColorModifier
from ovito.vis import Viewport, TachyonRenderer
import os

# 1. Load the trajectory (automatically detects LAMMPS dump format)
# Ensure "dump.lammpstrj" matches the output from your LAMMPS script
pipeline = import_file("dump.lammpstrj")

# Optional: Add a modifier to color atoms (e.g., standard blue)
pipeline.modifiers.append(AssignColorModifier(color=(0.1, 0.4, 0.8)))

# 2. Add pipeline to the visual scene
pipeline.add_to_scene()

# 3. Setup the viewport and camera
vp = Viewport()
vp.type = Viewport.Type.Perspective

# Automatically adjust the camera to fit all atoms in the simulation box
vp.zoom_all(size=(800, 600))

# 4. Determine the last frame and render
last_frame = pipeline.source.num_frames - 1

# Render using the high-quality Tachyon ray tracer
output_file = "lj_melt_render.png"
vp.render_image(
    filename=output_file,
    size=(1200, 900),
    frame=last_frame,
    background=(1.0, 1.0, 1.0), # White background
    renderer=TachyonRenderer(shadows=True, direct_light_intensity=1.2)
)

print(f"Successfully rendered frame {last_frame} to {output_file}")






# 1. Load the LAMMPS trajectory
pipeline = import_file("dump.lammpstrj")
pipeline.add_to_scene()

# 2. Configure the viewport and camera
vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.zoom_all(size=(800, 600))

# 3. Render the animation
# OVITO automatically iterates through all frames in the loaded trajectory
output_video = "lj_melt_animation.mp4"

print(f"Starting render of {pipeline.source.num_frames} frames...")

vp.render_anim(
    filename=output_video,
    size=(800, 600),
    fps=10, # Frames per second
    background=(1.0, 1.0, 1.0),
    renderer=TachyonRenderer(shadows=True, direct_light_intensity=1.2)
)

if os.path.exists(output_video):
    print(f"Render complete: {output_video}")
else:
    print("Render failed. Check if system video codecs are installed.")