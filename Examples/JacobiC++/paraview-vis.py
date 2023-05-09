###############################################################################
# Copyright (c) Lawrence Livermore National Security, LLC and other Ascent
# Project developers. See top-level LICENSE AND COPYRIGHT files for dates and
# other details. No copyright assignment is required to contribute to Ascent.
###############################################################################

# Same Python interpreter for all time steps
# We use count for one time initializations
try:
    count = count + 1
except NameError:
    count = 0

if count == 0:
    # Initialize ParaView
    import paraview
    paraview.options.batch = True
    paraview.options.symmetric = True
    from paraview.simple import LoadPlugin, Show, ColorBy, SelectPoints, \
        GetColorTransferFunction, GetActiveView, GetScalarBar, ResetCamera,\
        Render, SaveScreenshot, ExtractSelection, GetActiveCamera,\
        GetProperty, CreateRenderView, CreateWriter

    # load AscentSource to present simulation data as a VTK dataset
    LoadPlugin("/local/apps/Ascent/install/ascent-v0.9.1/examples/ascent/paraview-vis/paraview_ascent_source.py", remote=True, ns=globals())
    ascentSource = AscentSource()
    view = CreateRenderView()
    view.ViewSize = [1024, 1024]

    rep = Show()
    ColorBy(rep, ("POINTS", "temperature"))
    # rescale transfer function and show color bar
    transferFunction = GetColorTransferFunction('temperature')
    transferFunction.RescaleTransferFunction(0., 1.)
    renderView1 = GetActiveView()
    scalarBar = GetScalarBar(transferFunction, renderView1)
    scalarBar.Title = 'temperature'
    scalarBar.ComponentTitle = ''
    scalarBar.Visibility = 1
    rep.SetScalarBarVisibility(renderView1, True)

if not (count+1) % 1000:
  ascentSource.UpdateAscentData()
  ascentSource.UpdatePropertyInformation()
  cycle = GetProperty(ascentSource, "Cycle").GetElement(0)

  imageName = "image_{0:04d}.png".format(int(cycle))
  ResetCamera()
  Render()
  SaveScreenshot(imageName, ImageResolution=(1024, 1024))
  dataName = "jacobi_data_{0:04d}".format(int(cycle))
  writer = CreateWriter(dataName + ".pvti", ascentSource)
  writer.UpdatePipeline()


