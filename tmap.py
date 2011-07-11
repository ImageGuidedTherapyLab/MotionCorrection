import vtk
# echo vtk version info
print "using vtk version", vtk.vtkVersion.GetVTKVersion()
import vtk.util.numpy_support as vtkNumPy 
import numpy
import os
import scipy.io as scipyio

# raw data location
dirID =  115298
rootdir = "/FUS4/data2/CHUN_LI/070131/e114985/s%d" % dirID

dirID =  8980; nsteps = 60
rootdir = "/FUS4/data2/nanorods/20101019nanorods/e7605/s%d" % dirID

dirID =  75104; nsteps = 41
dirID =  75187; nsteps = 101
rootdir = "/FUS4/data2/2010_Gupta_Li_Rabbit/090612_rabbit2_9L033/e74953/s%d" % dirID

# get some header data
vtkDicomInfo = vtk.vtkDICOMImageReader()
vtkDicomInfo.SetFileName("%s/%s"%(rootdir,"i%d.MRDC.%d"%(dirID+1,1) ) )
vtkDicomInfo.Update()
dimensions = vtkDicomInfo.GetOutput().GetDimensions()
spacing_mm = vtkDicomInfo.GetOutput().GetSpacing()
origin_mm  = vtkDicomInfo.GetOutput().GetOrigin()
#convert to meter
spacing = [ 0.001 * dXi for dXi in spacing_mm  ]
origin  = [ 0.001 * Xi for Xi in origin_mm   ]
print dimensions, spacing, origin

# temperature map factor
alpha = +0.0097   # FIXME should be negative but phase messed up somewhere
# FIXME need to extract automagically at some point
#(0018|0081) Echo Time = 9.648
#(0018|0082) Inversion Time = 0
#(0018|0083) Number of Averages = 2
#(0018|0084) Imaging Frequency = 63.869849
echoTime = 9.648
echoTime = 7.84
imagFreq = 63.869849
tmap_factor = 1.0 / (2.0 * numpy.pi * imagFreq * alpha * echoTime * 1.e-3)
print tmap_factor

# return image data from raw file names
def GetRawDICOMData(filenames,fileID):
  print filenames,fileID
  vtkRealDcmReader = vtk.vtkDICOMImageReader()
  vtkRealDcmReader.SetFileName("%s/%s"%(rootdir,filenames[0]) )
  vtkRealDcmReader.Update()
  vtkRealData = vtk.vtkImageCast()
  vtkRealData.SetOutputScalarTypeToFloat()
  vtkRealData.SetInput( vtkRealDcmReader.GetOutput() )
  vtkRealData.Update( )
  real_image = vtkRealData.GetOutput().GetPointData() 
  real_array = vtkNumPy.vtk_to_numpy(real_image.GetArray(0)) 

  vtkImagDcmReader = vtk.vtkDICOMImageReader()
  vtkImagDcmReader.SetFileName("%s/%s"%(rootdir,filenames[1]) )
  vtkImagDcmReader.Update()
  vtkImagData = vtk.vtkImageCast()
  vtkImagData.SetOutputScalarTypeToFloat()
  vtkImagData.SetInput( vtkImagDcmReader.GetOutput() )
  vtkImagData.Update( )
  imag_image = vtkImagData.GetOutput().GetPointData() 
  imag_array = vtkNumPy.vtk_to_numpy(imag_image.GetArray(0)) 

  vtkAppend = vtk.vtkImageAppendComponents()
  vtkAppend.SetInput( 0,vtkRealDcmReader.GetOutput() )
  vtkAppend.SetInput( 1,vtkImagDcmReader.GetOutput() )
  vtkAppend.Update( )

  # write raw data
  vtkRawData = vtkAppend.GetOutput()
  vtkRawData.SetSpacing( spacing )
  vtkDcmWriter = vtk.vtkDataSetWriter()
  vtkDcmWriter.SetFileTypeToBinary()
  vtkDcmWriter.SetFileName("s%d/rawdata.%04d.vtk" % (dirID,fileID) )
  vtkDcmWriter.SetInput( vtkRawData )
  vtkDcmWriter.Update()

  # write raw phase data
  vtkPhase = vtk.vtkImageMathematics()
  vtkPhase.SetInput1( vtkRealData.GetOutput() )
  vtkPhase.SetInput2( vtkImagData.GetOutput() )
  vtkPhase.SetOperationToATAN2( )
  vtkPhase.Update( )
  vtkPhaseData = vtkPhase.GetOutput()
  vtkPhaseData.SetSpacing( spacing )
  vtkDcmWriter = vtk.vtkDataSetWriter()
  vtkDcmWriter.SetFileTypeToBinary()
  vtkDcmWriter.SetFileName("s%d/phase.%04d.vtk" % (dirID,fileID) )
  vtkDcmWriter.SetInput( vtkPhaseData)
  vtkDcmWriter.Update()

  return (real_array,imag_array)

# write a numpy data to disk in vtk format
def ConvertNumpyVTKImage(NumpyImageData):
  # Create initial image
  # imports raw data and stores it.
  dataImporter = vtk.vtkImageImport()
  # array is converted to a string of chars and imported.
  data_string = NumpyImageData.tostring()
  dataImporter.CopyImportVoidPointer(data_string, len(data_string))
  # The type of the newly imported data is set to unsigned char (uint8)
  dataImporter.SetDataScalarTypeToFloat()
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
  dataImporter.Update()
  return dataImporter.GetOutput()
  

# generate file names
realimagdata = []
for idfile in range(1,nsteps*2,2):
  realimagdata.append( ("i%d.MRDC.%d"%(dirID+idfile + 0,idfile + 0),
                        "i%d.MRDC.%d"%(dirID+idfile + 1,idfile + 1) )  ) 

# create working dir
os.system("mkdir -p s%d" % dirID)

deltat = 6.0
pvd=open("s%d/temperature.pvd" % dirID ,"w")
pvd.write('<?xml version="1.0"?>\n')
pvd.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
pvd.write('  <Collection>\n')
for idtime in range(nsteps):
     pvd.write('   <DataSet timestep="%f" part="0" file="%s.%04d.vti"/>\n' % (idtime*deltat,"temperature",idtime) )
pvd.write('  </Collection>\n')
pvd.write('</VTKFile>\n')

# create initial image as 1d array
absTemp = numpy.zeros(dimensions[0]*dimensions[1]*dimensions[2],
                       dtype=numpy.float32) + 21.0
vtkTempImage = ConvertNumpyVTKImage(absTemp)
vtkTempWriter = vtk.vtkXMLImageDataWriter()
vtkTempWriter.SetFileName( "s%d/temperature.%04d.vti" % (dirID,0))
vtkTempWriter.SetInput( vtkTempImage )
vtkTempWriter.Update()

# create a rendering window and renderer
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
 
# create a renderwindowinteractor
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
   
# loop and compute tmap
vtkPreviousImage = GetRawDICOMData( realimagdata.pop(0), 0  )
for idfile,dcmfiles in enumerate(realimagdata):
  
  # get current data set
  vtkCurrent_Image = GetRawDICOMData( dcmfiles , idfile + 1)

  #  - \delta \theta = atan( conj(S^i) * S^{i+1} ) 
  #                  = atan2(Im,Re) 
  #                  = atan2( S^{i+1}_y S^i_x - S^{i+1}_x S^i_y ,
  #                           S^{i+1}_x S^i_x + S^{i+1}_y S^i_y ) 
  deltaTemp = tmap_factor * numpy.arctan2(
                        vtkPreviousImage[0] * vtkCurrent_Image[1] 
                      - vtkPreviousImage[1] * vtkCurrent_Image[0] ,
                        vtkPreviousImage[1] * vtkCurrent_Image[1] 
                      + vtkPreviousImage[0] * vtkCurrent_Image[0]  )
  absTemp  = absTemp + deltaTemp 

  # write numpy to disk in vtk format
  vtkTempImage = ConvertNumpyVTKImage(absTemp)
  vtkTempWriter = vtk.vtkXMLImageDataWriter()
  vtkTempWriter.SetFileName( "s%d/temperature.%04d.vti" % (dirID,idfile+1))
  vtkTempWriter.SetInput( vtkTempImage )
  vtkTempWriter.Update()

  # write numpy to disk in matlab
  scipyio.savemat("s%d/temperature.%04d.mat" % (dirID,idfile+1), {'temp':absTemp})

  # update for next time step
  vtkPreviousImage = vtkCurrent_Image 

  # mapper
  mapper = vtk.vtkDataSetMapper()
  mapper.SetInput(vtkTempImage)
   
  # actor
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
   
  # assign actor to the renderer
  ren.AddActor(actor)
   
  # enable user interface interactor
  #iren.Initialize()
  #renWin.Render()
  #iren.Start()
