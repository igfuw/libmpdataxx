USAGE = 'Usage: pvpython plot.py timestep*.xmf'

from sys import argv, exit

try:
    from paraview.simple import *
except ImportError:
    print USAGE
    exit(1)

# checking arguments
if (
    len(argv) != 2 or
    not argv[1].startswith('timestep') or
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
view.CameraViewUp = [-0.26, -0.09, 0.96]
view.CameraFocalPoint = [61, 36, 32]
view.CameraViewAngle = 30
view.CameraPosition = [143, 69, 58]

# creating threshold of values > 0.5
th = Threshold()
th.ThresholdRange = [0.5, 4]

th_dp = GetDisplayProperties(th)
th_dp.Representation = "Surface With Edges"
th_dp.DiffuseColor = [1, 1, 1]
th_dp.EdgeColor = [0, 0, 0]

Show(th)

# domain diagonal
diag = Line()
diag.Point1 = [0, 0, 0]
diag.Point2 = [35, 35, 35]

diag_dp = GetDisplayProperties(diag)
diag_dp.AmbientColor = [0,0,0]
diag_dp.Representation = 'Wireframe'
diag_dp.LineWidth = 2
Show(diag)

# x axis
x_axis = Line()
x_axis.Point1 = [0, 0, 0]
x_axis.Point2 = [30, 0, 0]

x_axis_dp = GetDisplayProperties(x_axis)
x_axis_dp.AmbientColor = [0, 0, 0]
x_axis_dp.Representation = 'Wireframe'
x_axis_dp.LineWidth = 2
Show(x_axis)

# y axis
y_axis = Line()
y_axis.Point1 = [0, 0, 0]
y_axis.Point2 = [0, 30, 0]

y_axis_dp = GetDisplayProperties(y_axis)
y_axis_dp.AmbientColor = [0, 0, 0]
y_axis_dp.Representation = 'Wireframe'
y_axis_dp.LineWidth = 2
Show(y_axis)

# z axis
z_axis = Line()
z_axis.Point1 = [0, 0, 0]
z_axis.Point2 = [0, 0, 30]

z_axis_dp = GetDisplayProperties(z_axis)
z_axis_dp.AmbientColor = [0, 0, 0]
z_axis_dp.Representation = 'Wireframe'
z_axis_dp.LineWidth = 2
Show(z_axis)

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
