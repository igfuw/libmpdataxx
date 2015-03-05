USAGE = 'Usage: pvpython plot.py timestep*.xmf'

#TODO Paraview needs DISPLAY to be set, which won't work on Travis CI.
#We have a workaround in CMakeLists.txt that skips this test in such cases.
#It would be best to use the below code with cmake skip test but ...
# ... see https://public.kitware.com/Bug/view.php?id=13825
#import os, sys
#if (not os.environ.has_key("DISPLAY") or os.environ["DISPLAY"]==""):
#  sys.exit(44)

from sys import argv, exit
from math import sqrt

try:
    from paraview.simple import *
except ImportError:
    print USAGE
    exit(1)

# checking arguments
if (
    len(argv) != 2 or
    not argv[1].rpartition('/')[2].startswith('timestep') or
    not argv[1].endswith('.xmf')
   ):
    print USAGE
    exit(1)

# opening xmf file
fname = argv[1]
reader = OpenDataFile(fname)

# creating the view
view = CreateRenderView()
view.Background = [1, 1, 1]
view.CenterAxesVisibility = False
view.OrientationAxesVisibility = False
view.CameraViewUp = [0, 0, 1]

d = 25 / sqrt(3)
view.CameraFocalPoint = [50 - d, 50 + d, 50 + d]
view.CameraViewAngle = 45
view.CameraPosition = [200, 200, 100]

# creating threshold of values > 1
th = Threshold()
th.ThresholdRange = [1, 4]

th_dp = GetDisplayProperties(th)
th_dp.Representation = "Surface With Edges"
th_dp.DiffuseColor = [1, 1, 1]
th_dp.EdgeColor = [0, 0, 0]

Show(th)

L = 100

# z axis
line = Line()
line.Point1 = [0, 0, 0]
line.Point2 = [0, 0, L]
line_dp = GetDisplayProperties(line)
line_dp.AmbientColor = [0, 0, 0]
line_dp.Representation = 'Wireframe'
line_dp.LineWidth = 2
Show(line)

# xy plane
plane = Plane()
plane.Origin = [0, 0, 0]
plane.Point1 = [L, 0, 0]
plane.Point2 = [0, L, 0]
plane_dp = GetDisplayProperties(plane)
plane_dp.AmbientColor = [0, 0, 0]
plane_dp.Representation = 'Wireframe'
plane_dp.LineWidth = 2
Show(plane)

Render()

# export to SVG
exporters = servermanager.createModule('exporters')
if hasattr(exporters, 'GL2PSRenderViewExporterSVG'):
    svg_exporter = exporters.GL2PSRenderViewExporterSVG()
else:
    svg_exporter = exporters.GL2PSExporterSVG()
svg_exporter.Rasterize3Dgeometry = False;
svg_exporter.GL2PSdepthsortmethod = 'BSP sorting (slow, best)';
svg_exporter.SetView(view)
svg_exporter.FileName = fname.rpartition('.')[0] + ".svg"
svg_exporter.Write()
