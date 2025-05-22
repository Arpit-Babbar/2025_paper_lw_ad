# state file generated using paraview version 5.10.1

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [500, 500]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 2.7755575615628914e-17, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.3984907353079417, 0.025336687097089222, 10000.0]
renderView1.CameraFocalPoint = [0.3984907353079417, 0.025336687097089222, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.39904972451904025

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(500, 500)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'XML Rectilinear Grid Reader'
import os
dir_path = os.path.dirname(os.path.realpath(__file__))
filename = os.path.join(dir_path, 'output_m2000', 'sol100.vtr')
sol_vtr = XMLRectilinearGridReader(registrationName='sol100.vtr', FileName=[filename])
sol_vtr.PointArrayStatus = ['sol']
sol_vtr.TimeArray = 'None'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sol100vtr
sol100vtrDisplay = Show(sol_vtr, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'sol'
solLUT = GetColorTransferFunction('sol')
solLUT.RGBPoints = [0.015462454269614007, 0.0, 0.0, 0.34902, 0.01962313400107963, 0.039216, 0.062745, 0.380392, 0.02490338088055279, 0.062745, 0.117647, 0.411765, 0.03160445111610415, 0.090196, 0.184314, 0.45098, 0.04010866376501589, 0.12549, 0.262745, 0.501961, 0.0509012133482482, 0.160784, 0.337255, 0.541176, 0.06459785186321214, 0.2, 0.396078, 0.568627, 0.08198001954869137, 0.239216, 0.454902, 0.6, 0.1040394287326327, 0.286275, 0.521569, 0.65098, 0.13203464442434792, 0.337255, 0.592157, 0.701961, 0.16756288976811692, 0.388235, 0.654902, 0.74902, 0.2126511731057798, 0.466667, 0.737255, 0.819608, 0.26987193576001867, 0.572549, 0.819608, 0.878431, 0.3424898186413068, 0.654902, 0.866667, 0.909804, 0.43464792121720475, 0.752941, 0.917647, 0.941176, 0.5516041795575074, 0.823529, 0.956863, 0.968627, 0.7000313496340426, 0.988235, 0.960784, 0.901961, 0.7000313496340426, 0.941176, 0.984314, 0.988235, 0.8153630977943276, 0.988235, 0.945098, 0.85098, 0.9496960123176063, 0.980392, 0.898039, 0.784314, 1.1274502109618318, 0.968627, 0.835294, 0.698039, 1.4308276152982153, 0.94902, 0.733333, 0.588235, 1.8158386461726306, 0.929412, 0.65098, 0.509804, 2.3044495043847926, 0.909804, 0.564706, 0.435294, 2.92453711647376, 0.878431, 0.458824, 0.352941, 3.7114796090600324, 0.839216, 0.388235, 0.286275, 4.710174752398977, 0.760784, 0.294118, 0.211765, 5.97760153227826, 0.701961, 0.211765, 0.168627, 7.586070996728218, 0.65098, 0.156863, 0.129412, 9.627351849508017, 0.6, 0.094118, 0.094118, 12.217906169636386, 0.54902, 0.066667, 0.098039, 15.505533972737005, 0.501961, 0.05098, 0.12549, 19.67780570922952, 0.45098, 0.054902, 0.172549, 24.97276380233136, 0.4, 0.054902, 0.192157, 31.69250378534457, 0.34902, 0.070588, 0.211765]
solLUT.UseLogScale = 1
solLUT.ColorSpace = 'Lab'
solLUT.NanColor = [0.25, 0.0, 0.0]
solLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'sol'
solPWF = GetOpacityTransferFunction('sol')
solPWF.Points = [0.015462454269614003, 0.0, 0.5, 0.0, 31.692503785344584, 1.0, 0.5, 0.0]
solPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
sol100vtrDisplay.Representation = 'Surface'
sol100vtrDisplay.ColorArrayName = ['POINTS', 'sol']
sol100vtrDisplay.LookupTable = solLUT
sol100vtrDisplay.SelectTCoordArray = 'None'
sol100vtrDisplay.SelectNormalArray = 'None'
sol100vtrDisplay.SelectTangentArray = 'None'
sol100vtrDisplay.OSPRayScaleArray = 'sol'
sol100vtrDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
sol100vtrDisplay.SelectOrientationVectors = 'None'
sol100vtrDisplay.ScaleFactor = 0.09997654496148467
sol100vtrDisplay.SelectScaleArray = 'None'
sol100vtrDisplay.GlyphType = 'Arrow'
sol100vtrDisplay.GlyphTableIndexArray = 'None'
sol100vtrDisplay.GaussianRadius = 0.004998827248074234
sol100vtrDisplay.SetScaleArray = ['POINTS', 'sol']
sol100vtrDisplay.ScaleTransferFunction = 'PiecewiseFunction'
sol100vtrDisplay.OpacityArray = ['POINTS', 'sol']
sol100vtrDisplay.OpacityTransferFunction = 'PiecewiseFunction'
sol100vtrDisplay.DataAxesGrid = 'GridAxesRepresentation'
sol100vtrDisplay.PolarAxes = 'PolarAxesRepresentation'
sol100vtrDisplay.ScalarOpacityUnitDistance = 0.008909867778506182
sol100vtrDisplay.ScalarOpacityFunction = solPWF
sol100vtrDisplay.OpacityArrayName = ['POINTS', 'sol']
sol100vtrDisplay.SliceFunction = 'Plane'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sol100vtrDisplay.ScaleTransferFunction.Points = [0.015462454269614003, 0.0, 0.5, 0.0, 31.692503785344584, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sol100vtrDisplay.OpacityTransferFunction.Points = [0.015462454269614003, 0.0, 0.5, 0.0, 31.692503785344584, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
sol100vtrDisplay.SliceFunction.Origin = [0.5, 2.7755575615628914e-17, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for solLUT in view renderView1
solLUTColorBar = GetScalarBar(solLUT, renderView1)
solLUTColorBar.Orientation = 'Horizontal'
solLUTColorBar.WindowLocation = 'Any Location'
solLUTColorBar.Position = [0.07199999999999983, 0.8780000000000002]
solLUTColorBar.Title = ''
solLUTColorBar.ComponentTitle = ''
solLUTColorBar.TitleFontSize = 20
solLUTColorBar.LabelFontSize = 20
solLUTColorBar.RangeLabelFormat = '%2.2f'
solLUTColorBar.ScalarBarLength = 0.8240000000000003

# set color bar visibility
solLUTColorBar.Visibility = 1

# show color legend
sol100vtrDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# restore active source
SetActiveSource(sol_vtr)
# ----------------------------------------------------------------

myview = GetActiveView()
out_filename = os.path.join(dir_path, 'paper_figures', 'm2000.png')
SaveScreenshot(out_filename, myview,
        ImageResolution=[1500, 1500])


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
