"""
background correction model for 2d imaging
    IMPORTANT Points: 
      -  Numpy stores arrays in c-style row-major order,
         VTK/Matlab is expected fortran-style column major order
         the order='F' flag is needed to convert between the styles   
      -  The = operator for numpy arrays just copies by reference
         .copy() provides a deep copy (ie physical memory copy)
"""

# import needed modules
import petsc4py, numpy, sys,os
import scipy.io as scipyio
PetscOptions = sys.argv
PetscOptions.append("-ksp_monitor") 
PetscOptions.append("-ksp_rtol") 
PetscOptions.append("1.e-12") 
# need to solve w/ constant null space for Neumann problem
#PetscOptions.append("-ksp_constant_null_space")
#PetscOptions.append("-help")
petsc4py.init(PetscOptions)

# get rank
from petsc4py import PETSc
petscRank = PETSc.COMM_WORLD.getRank()
petscSize = PETSc.COMM_WORLD.Get_size()
sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

# set shell context
# TODO import vtk should be called after femLibrary ???? 
# FIXME WHY IS THIS????
import femLibrary
# initialize libMesh data structures
libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  
# load vtk modules to read imaging
import vtk 
import vtk.util.numpy_support as vtkNumPy 

#FileNameTemplate = "/data/fuentes/biotex/MayoLiverDataJan2011/VTKData/226p-536/phase.%04d.vti"
#FileNameTemplate = "/data/fuentes/mdacc/090612_rabbit2_9L033/s75187/phase.%04d.vtk"
#FileNameTemplate = "/home/jyung/e137/Processed/s22000/phase.%04d.vtk"
FileNameTemplate = "/data/jyung/old/motion_phantom_july11/multiplanar/skip2/phaseROI.%04d.vtk"
 
# set the default reader based on extension
if( FileNameTemplate.split(".").pop() == "vtk"):
   vtkImageReader = vtk.vtkDataSetReader
elif( FileNameTemplate.split(".").pop() == "vti"):
   vtkImageReader = vtk.vtkXMLImageDataReader
else:
   raise RuntimeError("uknown file")

# get dimension info from header
vtkSetupReader = vtkImageReader() 
vtkSetupReader.SetFileName(FileNameTemplate % 0 ) 
vtkSetupReader.Update() 
dimensions = vtkSetupReader.GetOutput().GetDimensions()
numberPointsImage =  vtkSetupReader.GetOutput().GetNumberOfPoints()
spacing = vtkSetupReader.GetOutput().GetSpacing()
origin  = vtkSetupReader.GetOutput().GetOrigin()
print "#points",numberPointsImage , "dimensions ",dimensions , "spacing ",spacing , "origin ",origin 
# temperature map factor
alpha = +0.0097    # FIXME should be neg
# FIXME need to extract automagically at some point
#(0018|0081) Echo Time = 9.648
#(0018|0082) Inversion Time = 0
#(0018|0083) Number of Averages = 2
#(0018|0084) Imaging Frequency = 63.869849
echoTime = 9.648
imagFreq = 63.869849
tmap_factor = - 1.0 / (2.0 * numpy.pi * imagFreq * alpha * echoTime * 1.e-3)
print "tmap_factor = ", tmap_factor

# return image data from raw file names
def GetNumpyPhaseData(filename):
  vtkReader = vtkImageReader() 
  vtkReader.SetFileName(filename) 
  vtkReader.Update() 
  vtkImageCast = vtk.vtkImageCast()
  vtkImageCast.SetOutputScalarTypeToDouble()
  vtkImageCast.SetInput( vtkReader.GetOutput() )
  vtkImageCast.Update( )
  phase_image = vtkImageCast.GetOutput().GetPointData() 
  # convert to phase
  phase_array =  vtkNumPy.vtk_to_numpy(phase_image.GetArray(0)) 
  #return (2.0 * numpy.pi / 4095. ) *phase_array 
  return phase_array.reshape(dimensions,order='F')

# write a numpy data to disk in vtk format
def ConvertNumpyVTKImage(NumpyImageData, arrayName  ):
  # Create initial image
  # imports raw data and stores it.
  dataImporter = vtk.vtkImageImport()
  # array is converted to a string of chars and imported.
  # numpy array stored as ROW MAJOR
  # MUST write out in COLUMN MAJOR format to be the same as VTK
  data_string = NumpyImageData.tostring(order='F')
  dataImporter.CopyImportVoidPointer(data_string, len(data_string))
  # The type of the newly imported data is set to unsigned char (uint8)
  dataImporter.SetDataScalarTypeToDouble()
  # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
  # must be told this is the case.
  dataImporter.SetNumberOfScalarComponents(1)
  # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
  # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
  # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
  # VTK complains if not both are used.
  dataImporter.SetDataExtent( 0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  dataImporter.SetWholeExtent(0, dimensions[0]-1, 0, dimensions[1]-1, 0, dimensions[2]-1)
  dataImporter.SetDataSpacing( spacing )
  dataImporter.SetDataOrigin(  origin  )
  dataImporter.SetScalarArrayName(  arrayName  )
  dataImporter.Update()
  return dataImporter.GetOutput()

# store control variables
getpot = femLibrary.PylibMeshGetPot(PetscOptions) 

# initialize FEM Mesh
femMesh = femLibrary.PylibMeshMesh()
#ROI = [[100,150],   # pixel # of xbounds
#       [100,150],   # pixel # of ybounds
#       [  0,29]]   # pixel # of zbounds
#ROI = [[50,200],   # pixel # of xbounds
#       [50,200],   # pixel # of ybounds
#       [ 0,29]]   # pixel # of zbounds
ROI = [[0,150],[50,200],[0,29]]
npixelROI = tuple( [ (pixel[1] - pixel[0] ) for pixel in ROI] )
nelemROI  = [ (pixel[1] - pixel[0] - 1 ) for pixel in ROI] 
if( nelemROI[2] < 0  ):
   nelemROI[2] = 0 
nelemROI  = tuple( nelemROI )
xbounds = [ origin[0]+spacing[0]*(ROI[0][0]+0.5),origin[0]+spacing[0]*(ROI[0][1]+0.5) ]
ybounds = [ origin[1]+spacing[1]*(ROI[1][0]+0.5),origin[1]+spacing[1]*(ROI[1][1]+0.5) ]
zbounds = [ origin[2]+spacing[2]*(ROI[2][0]+0.5),origin[2]+spacing[2]*(ROI[2][1]+0.5) ]
# setup structure Grid expect # of elements
femMesh.SetupStructuredGrid( nelemROI , 
                             xbounds,ybounds,zbounds,[1,1,1,1,1,1]) 

# setup imaging to interpolate onto FEM mesh
femImaging = femLibrary.PytttkImaging(getpot, dimensions ,origin,spacing) 

# add the data structures for the Background System Solve
eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
bgSystem   = eqnSystems.AddBackgroundSystem( "Background" ) 
maskSystem = eqnSystems.AddExplicitSystem( "ImageMask" ) 
maskSystem.AddFirstLagrangeVariable( "mask" ) 
# initialize libMesh data structures
eqnSystems.init( ) 
# print info
eqnSystems.PrintSelf() 
  
# create working dir
os.system("mkdir -p Processed")

# write IC
MeshOutputFile = "Processed/fem_data.e"
exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  

# loop over desired time instances
ntime=0
for timeID in range(0,ntime+1):
   print "working on time id %d " % timeID
   # read in data
   phase_curr = GetNumpyPhaseData(FileNameTemplate % timeID ) 
   
   # write original delta tmap 
   vtkTempImage = ConvertNumpyVTKImage(phase_curr,"orthogphase")
   vtkTmpWriter = vtk.vtkDataSetWriter()
   vtkTmpWriter.SetFileName("Processed/orthogphase.%04d.vtk" % timeID )
   vtkTmpWriter.SetInput(vtkTempImage)
   vtkTmpWriter.Update()

   vtkImageMask = None
   #default look for a user defined image mask
   maskFileName = "%s/SNRuncert.%04d.vtk" %  (
                 FileNameTemplate[:FileNameTemplate.rfind("/")],timeID)
   if( os.path.isfile( maskFileName ) ):
     vtkMaskReader = vtk.vtkDataSetReader() 
     vtkMaskReader.SetFileName( maskFileName ) 
     vtkMaskReader.Update() 
     # get image mask
     vtkImageThresh = vtk.vtkImageThreshold() 
     vtkImageThresh.ReplaceOutOn()
     vtkImageThresh.SetOutValue(0.0)
     vtkImageThresh.ReplaceInOn()
     vtkImageThresh.SetInValue(1.0)
     vtkImageThresh.SetOutputScalarTypeToDouble()
     # set threshold
     vtkImageThresh.ThresholdByLower( 100.0)
     vtkImageThresh.SetInput(vtkMaskReader.GetOutput())
     #else: # use gradient as threshold
     #  vtkPhaseData = ConvertNumpyVTKImage( phase_curr )
     #  vtkGradientImage = vtk.vtkImageGradient() 
     #  vtkGradientImage.SetInput(vtkPhaseData) 
     #  vtkGradientImage.Update() 
     #  
     #  # take magnitude of gradient
     #  vtkImageNorm = vtk.vtkImageMagnitude() 
     #  vtkImageNorm.SetInput(vtkGradientImage.GetOutput())
     #  vtkImageNorm.Update() 
     #  # set threshold
     #  vtkImageThresh.ThresholdByLower( 100.0* tmap_factor * 2.0*numpy.pi/4095.)
     #  vtkImageThresh.SetInput(vtkImageNorm.GetOutput())
     vtkImageThresh.Update( )
     vtkImageMask = vtkImageThresh.GetOutput( )
   else: # if nothing available use 2D ROI as a mask
     # vtk uses fortran style storage
     numpyImageMask = 1.0 + numpy.zeros(dimensions,
                                        dtype=numpy.double,order='F')
     # TODO: numpy seems to index in the transpose of VTK ? 
     numpyImageMask[(ROI[0][0]+1):(ROI[0][1]-1),
                    (ROI[1][0]+1):(ROI[1][1]-1),:] = 0.0
     numpyImageMask.tofile("Processed/threshold.%04d.dat" % timeID )
     vtkImageMask = ConvertNumpyVTKImage( numpyImageMask , "theshold")
   # check output
   vtkWriter = vtk.vtkXMLImageDataWriter()
   vtkWriter.SetFileName("Processed/threshold.%04d.vti" % timeID )
   vtkWriter.SetInput( vtkImageMask )
   vtkWriter.Update()

   # pass pointer to c++
   image_cells = vtkImageMask.GetPointData() 
   data_array = vtkNumPy.vtk_to_numpy( image_cells.GetArray(0) ) 
   # need to pass numpy array's w/ Fortran storage... ie painful to debug
   v1 = PETSc.Vec().createWithArray(numpy.ravel(phase_curr,order='F'), comm=PETSc.COMM_SELF)
   v2 = PETSc.Vec().createWithArray(data_array, comm=PETSc.COMM_SELF)

   # Project imaging onto libMesh data structures
   femImaging.ProjectImagingToFEMMesh("Background",0.0,v1,eqnSystems)  
   femImaging.ProjectImagingToFEMMesh("ImageMask" ,0.0,v2,eqnSystems)  
   bgSystem.SystemSolve( ) 
   exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, timeID+1, timeID )  

   # get libMesh Background Solution as numpy data structure
   maxwell_array = eqnSystems.GetSolutionVector( "Background" )[...]
   maxwell_data  = phase_curr.copy()
   # reshape from colume major Fortran-like storage
   maxwell_data[ROI[0][0]:ROI[0][1],
                ROI[1][0]:ROI[1][1],
                ROI[2][0]:ROI[2][1]] = maxwell_array.reshape( npixelROI ,order='F')
   # write numpy to disk in matlab
   scipyio.savemat("Processed/background.%04d.mat"%(timeID), {'maxwell':maxwell_array} )
   # check output
   vtkTempImage = ConvertNumpyVTKImage(maxwell_data,"maxwell")
   vtkWriterTmpTwo = vtk.vtkXMLImageDataWriter()
   vtkWriterTmpTwo.SetFileName("Processed/maxwell.%04d.vti" % timeID )
   vtkWriterTmpTwo.SetInput( vtkTempImage )
   vtkWriterTmpTwo.Update()

   # compute temperature difference
   delta_temp =tmap_factor * ( phase_curr - maxwell_data ) 

   # check output
   vtkTempImage = ConvertNumpyVTKImage(delta_temp,"deltat")
   vtkWriterTmpTwo = vtk.vtkXMLImageDataWriter()
   vtkWriterTmpTwo.SetFileName("Processed/deltat.%04d.vti" % timeID )
   vtkWriterTmpTwo.SetInput( vtkTempImage )
   vtkWriterTmpTwo.Update()

