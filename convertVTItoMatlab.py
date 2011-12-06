
####################################################################
def ConvertVTItoMatlab(work_dir):
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  import numpy
  import scipy.io as scipyio
  import os
  # echo vtk version info
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()

  directoryList = os.listdir(work_dir)
  files = filter(lambda x:os.path.isfile("%s/%s" % (work_dir,x) ),directoryList)
  vtiFiles = filter(lambda x: x.split(".").pop() == "vti" ,files)
  filenames = filter(lambda x: x.split(".").pop(0) == "updatefemROI5" ,vtiFiles)
  
  for idlist, fileID in enumerate( filenames):
    print "reading %s/%s" % (work_dir,fileID) 
    vtkReader = vtk.vtkXMLImageDataReader() 
    vtkReader.SetFileName( "%s/%s" % (work_dir,fileID) ) 
    vtkReader.Update()
    imageDataVTK = vtkReader.GetOutput()
    dimensions = imageDataVTK.GetDimensions()
    spacing = imageDataVTK.GetSpacing()
    origin  = imageDataVTK.GetOrigin()
    print spacing, origin, dimensions
    #fem.SetImagingDimensions( dimensions ,origin,spacing) 

    image_point_data = imageDataVTK.GetPointData() 
    image_data       = vtkNumPy.vtk_to_numpy( image_point_data.GetArray(0) ) 
    # write numpy to disk in matlab
    scipyio.savemat("%s/%s.mat" % (work_dir,fileID), {'spacing':spacing, 'origin':origin,'image':image_data})
    print "wrote %d of %d" % (idlist, len( filenames))

  
# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--work_dir",
                  action="store", dest="work_dir", default=None,
                  help="project data in exodus FILE  DIR/fem_stats.e to imaging", metavar = "DIR")
(options, args) = parser.parse_args()
if (options.work_dir):
  ConvertVTItoMatlab(options.work_dir)
else:
  parser.print_help()
  print options
