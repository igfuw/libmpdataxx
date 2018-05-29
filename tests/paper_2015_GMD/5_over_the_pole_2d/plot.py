USAGE = 'Usage: pvpython plot.py timestep*.xmf'

from sys import argv, exit

try:
    from paraview.simple import *
except ImportError:
    print USAGE
    exit(1)

from math import pi
from numpy import linspace
import os

# checking arguments
if (
    len(argv) != 2 or
    #not argv[1].startswith('timestep') or
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

# ortographic projection
ort = Calculator()
ort.CoordinateResults = True

dx = pi / 64
dy = pi / 64

x0 = 3 * pi / 2
y0 = 0

sx = '(coordsX - 0.5) * {dx}'.format(dx = dx)
sy = '(coordsY) * {dy} - {pi} / 2'.format(dy = dy, pi = pi)

ort.Function = '''
                 cos({y}) * sin({x} - {x0}) * iHat
               + (cos({y0}) * sin({y}) - sin({y0}) * cos({y}) * cos({x} - {x0})) * jHat
               + kHat
               '''.format(x = sx, x0 = x0, y = sy, y0 = y0)

ort_dp = GetDisplayProperties(ort)
ort_dp.Representation = "Surface"
ort_dp.ColorArrayName = 'psi'

# assigning color map to scalar values
psi = ort.CellData.GetArray('psi')
path = os.path.abspath(os.path.dirname(__file__))

if "LoadLookupTable" in dir(paraview.simple):
	LoadLookupTable(path + '/cmap.xml')
	ort_dp.LookupTable=AssignLookupTable(psi, 'my_set1', [0, 1])

# psi contours
ct = Contour()
ct.ContourBy = 'psi'
ct.Isosurfaces = [0.1 * i for i in xrange(10)]

ct_dp = GetDisplayProperties(ct)
ct_dp.Representation = "Wireframe"
ct_dp.AmbientColor = [0, 0, 0]

Show(ct)

# x coordinate lines
SetActiveSource(ort)
x_cl = Calculator()

rho = 'sqrt(coordsX * coordsX + coordsY * coordsY)'
c = 'asin({rho})'.format(rho = rho)

x_cl.Function = '{x0} + atan(coordsX * sin({c}) / ({rho} * cos({c})))'.format(x0\
 = x0, rho = rho, c = c)
x_cl.ResultArrayName = 'x'

x_ct = Contour()
x_ct.ContourBy = 'x'
iso = [3.17, 6.24]
iso.extend(linspace(3.17, 6.24, 8))
x_ct.Isosurfaces = iso

x_ct_dp = GetDisplayProperties(x_ct)
x_ct_dp.Representation = "Wireframe"
x_ct_dp.AmbientColor = [0, 0, 0]
x_ct_dp.LineWidth = 1

Show(x_ct)

# y coordinate lines
SetActiveSource(ort)

y_cl = Calculator()
y_cl.Function = 'asin(coordsY * sin({c}) / {rho})'.format(rho = rho, c = c)
y_cl.ResultArrayName = 'y'

y_ct = Contour()
y_ct.ContourBy = 'y'
iso = []
iso.extend(linspace(-pi / 2, pi / 2, 9))
y_ct.Isosurfaces = iso

y_ct_dp = GetDisplayProperties(y_ct)
y_ct_dp.Representation = "Wireframe"
y_ct_dp.AmbientColor = [0, 0, 0]
y_ct_dp.LineWidth = 1

# color bar
if "LoadLookupTable" in dir(paraview.simple):
	lt = GetLookupTableForArray('psi', 6)
	bar = CreateScalarBar(LookupTable = ort_dp.LookupTable)
	bar.Position = [0.9, 0.25]
	bar.LabelColor = [0, 0, 0]
	bar.NumberOfLabels = 6
	bar.LabelFontSize = 16
	view.Representations.append(bar)

Render()

# svg export
exporters = servermanager.createModule('exporters')
if hasattr(exporters, 'GL2PSRenderViewExporterSVG'):
    svg_exporter = exporters.GL2PSRenderViewExporterSVG()
else:
    svg_exporter = exporters.GL2PSExporterSVG()
svg_exporter.GL2PSdepthsortmethod = 'BSP sorting (slow, best)';
svg_exporter.Rasterize3Dgeometry = False;
svg_exporter.SetView(view)
svg_exporter.FileName = fname.rpartition('.')[0] + ".svg"
svg_exporter.Write()
