CMake_Minimum_Required(VERSION 2.8)

Project(Deformetrica)
Set (DEFORMETRICA_VERSION_MAJOR 2)
Set (DEFORMETRICA_VERSION_MINOR 1)

If (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)


enable_testing()

#
# Find packages (Armadillo, CUDA, ITK & VTK) :
#
# Armadillo :
# https://cmake.org/cmake/help/v3.0/module/FindArmadillo.html
Find_Package(Armadillo REQUIRED)
Include_Directories(${ARMADILLO_INCLUDE_DIRS})
IF (ARMADILLO_FOUND)
    Include_directories(${ARMADILLO_INCLUDE_DIRS})
    Message("Armadillo is found")
ELSE (ARMADILLO_FOUND)
    MESSAGE(FATAL_ERROR "Armadillo was not found :(")
ENDIF (ARMADILLO_FOUND)


# ITK :
Find_Package(ITK)
If (ITK_FOUND)
    Include (${ITK_USE_FILE})
Else (ITK_FOUND)
    Message (FATAL_ERROR "ITK not found. Please set ITK_DIR.")
Endif (ITK_FOUND)

# VTK :
Find_Package(VTK)
If (VTK_FOUND)
    If ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
        Include (${USE_VTK_FILE})
    Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Find_Package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkWrappingTools vtkIOLegacy vtkFiltersCore vtkFiltersGeneral NO_MODULE)
        If (VTK_FOUND)
            Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
            Include (${VTK_USE_FILE})
        Else (VTK_FOUND)
            Message (FATAL_ERROR "Some components of VTK have not been found.")
        Endif (VTK_FOUND)
    Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )
Else (VTK_FOUND)
    Message (FATAL_ERROR "VTK not found. Please set VTK_DIR.")
Endif (VTK_FOUND)

# Cuda :
Find_Package(CUDA)
If (CUDA_FOUND)
    Option (USE_CUDA "Use cuda implementation of convolutions" ON) 
Else (CUDA_FOUND)
    Message ("CUDA was not found.")
Endif (CUDA_FOUND)

# CMRep
SET( CMREP_SOURCE_DIR ${CMREP_SOURCE_DIR} ../cmrep )

MESSAGE( ${CMREP_SOURCE_DIR} )

# Pardiso : CMRep
# Optionally, support PARDISO
OPTION(USE_PARDISO "Use the PARDISO sparse solver library" ON)
IF(USE_PARDISO)
    # Search for PARDISO
    FIND_LIBRARY(PARDISO_LIB pardiso DOC "The location of the PARDISO shared library")
    IF(PARDISO_LIB)
      ADD_DEFINITIONS(-DHAVE_PARDISO)
      SET(SOLVER_SRC ${CMREP_SOURCE_DIR}/src/PardisoInterface.cxx)
    ENDIF(PARDISO_LIB)  

    # On Linux, Lapack is also required
    INCLUDE(${CMREP_SOURCE_DIR}/CMake/FindLapack.cmake)

    FIND_LIBRARY(GOMP_LIB gomp DOC "The location of the libgomp library")

    # Set the PDE libraries variable
    IF(PARDISO_LIB AND LAPACK_FOUND)
      SET(PDE_LIBS ${PARDISO_LIB} ${LAPACK_LIBS} ${GOMP_LIB})
    ENDIF(PARDISO_LIB AND LAPACK_FOUND)
ELSE(USE_PARDISO)
  ADD_DEFINITIONS(-UHAVE_PARDISO)
ENDIF(USE_PARDISO)

# Set CMRep include directories
INCLUDE_DIRECTORIES(${CMREP_SOURCE_DIR}/src ${CMREP_SOURCE_DIR}/src/dijkstra)

# Get GLUT
FIND_PACKAGE(GLUT)
IF(GLUT_FOUND)
  INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
  LINK_DIRECTORIES(${GLUT_INCLUDE_DIR})
ENDIF(GLUT_FOUND)

FIND_PACKAGE(OpenGL)
IF(OpenGL_FOUND)

ENDIF(OpenGL_FOUND)

# Include subdirectories
SUBDIRS(../cmrep/extras/toms611)

# Sources for the PDE executable
SET(COMMON_SRCS
  ${CMREP_SOURCE_DIR}/src/BasisFunctions2D.cxx
  ${CMREP_SOURCE_DIR}/src/BranchingSubdivisionSurface.cxx
  ${CMREP_SOURCE_DIR}/src/BruteForceSubdivisionMedialModel.cxx
  ${CMREP_SOURCE_DIR}/src/CoefficientMapping.cxx
  ${CMREP_SOURCE_DIR}/src/CartesianMedialModel.cxx
  ${CMREP_SOURCE_DIR}/src/DiffeomorphicEnergyTerm.cxx
  ${CMREP_SOURCE_DIR}/src/GeometryDescriptor.cxx
  ${CMREP_SOURCE_DIR}/src/ITKImageWrapper.cxx
  ${CMREP_SOURCE_DIR}/src/JacobianDistortionPenaltyTerm.cxx
  ${CMREP_SOURCE_DIR}/src/MedialAtom.cxx
  ${CMREP_SOURCE_DIR}/src/MedialAtomGrid.cxx
  ${CMREP_SOURCE_DIR}/src/MeshMedialPDESolver.cxx
  ${CMREP_SOURCE_DIR}/src/MedialModelIO.cxx
  ${CMREP_SOURCE_DIR}/src/MedialPDEMasks.cxx
  ${CMREP_SOURCE_DIR}/src/MedialPDESites.cxx
  ${CMREP_SOURCE_DIR}/src/MeshTraversal.cxx
  ${CMREP_SOURCE_DIR}/src/OptimizationTerms.cxx
  ${CMREP_SOURCE_DIR}/src/OptimizationParameters.cxx
  ${CMREP_SOURCE_DIR}/src/PDESubdivisionMedialModel.cxx
  ${CMREP_SOURCE_DIR}/src/PrincipalComponents.cxx
  ${CMREP_SOURCE_DIR}/src/PrincipalComponentsPenaltyTerm.cxx
  ${CMREP_SOURCE_DIR}/src/Procrustes.cxx
  ${CMREP_SOURCE_DIR}/src/Registry.cxx
  ${CMREP_SOURCE_DIR}/src/ScriptImaging.cxx
  ${CMREP_SOURCE_DIR}/src/ScriptInterface.cxx
  ${CMREP_SOURCE_DIR}/src/SmoothedImageSampler.cxx
  ${CMREP_SOURCE_DIR}/src/SparseMatrix.cxx
  ${CMREP_SOURCE_DIR}/src/SparseSolver.cxx
  ${CMREP_SOURCE_DIR}/src/SubdivisionMedialModel.cxx
  ${CMREP_SOURCE_DIR}/src/SubdivisionSurface.cxx
  ${CMREP_SOURCE_DIR}/src/System.cxx
  ${CMREP_SOURCE_DIR}/src/TestSolver.cxx
  ${CMREP_SOURCE_DIR}/src/VTKMedialMesh.cxx
  ${CMREP_SOURCE_DIR}/src/VTKSubdivision.cxx
  ${SOLVER_SRC})

# Create the CM-REP library
ADD_LIBRARY(Deformetrica ${COMMON_SRCS})

# Deformetrica
Include_Directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/atlases
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/4D_atlases
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/deformable_objects
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/deformations
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/energy_minimization
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/kernels
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/linear_algebra
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/algo/utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/io
    ${PROJECT_BINARY_DIR}
)

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
    Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)

# Basic flags for the C++ compiler :
#Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -llapack")
#Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -larmadillo")

# Extra flags for the C++ compiler :
If (CMAKE_BUILD_TYPE MATCHES Release)
    Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DARMA_NO_DEBUG -DNDEBUG -msse4.2")
Endif (CMAKE_BUILD_TYPE MATCHES Release)

# Extra flags for the C++ compiler if version is Gprof :
If (CMAKE_BUILD_TYPE MATCHES Gprof)
    Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
Endif (CMAKE_BUILD_TYPE MATCHES Gprof)

# Set computations on double precision by default :
If (NOT DEFINED SET_DOUBLE_PRECISION)
    Set (SET_DOUBLE_PRECISION YES)
Endif (NOT DEFINED SET_DOUBLE_PRECISION)
If (SET_DOUBLE_PRECISION)
    Option (USE_DOUBLE_PRECISION "Deformetrica will make computations in double precision" ON) 
Else (SET_DOUBLE_PRECISION)
    Message ("Computations will be in simple precision")
Endif (SET_DOUBLE_PRECISION)


# Configure a header file to pass some of the CMake settings to the source code :
Configure_File (
    "${PROJECT_SOURCE_DIR}/DeformetricaConfig.h.in"
    "${PROJECT_BINARY_DIR}/DeformetricaConfig.h"
)

If (USE_CUDA)
    Message("Deformetrica will be compiled with Cuda")
else(USE_CUDA)
    Message("Deformetrica will be compiled without Cuda")
    Add_Executable(sparseGeodesicRegression3CMRep
	  ../src/io/SparseDiffeoParametersXMLFile.cxx
	  ../src/io/DeformableObjectParametersXMLFile.cxx
	  ../src/io/SparseDiffeoParameters.cxx
	  ../src/io/DeformableObjectParameters.cxx
	  ../src/algo/utilities/SimpleTimer.cxx
	  ../src/algo/utilities/myvtkPolyDataNormals.cxx
	  ${COMMON_SRCS}
	  regression/sparseRegression3CMRep.cxx)
endif(USE_CUDA)


# Add_Executable(CreateDeformationField
#   ../src/io/SparseDiffeoParametersXMLFile.cxx
#   ../src/io/DeformationFieldIO.txx
#   ../src/io/SparseDiffeoParameters.cxx
#   ../src/algo/utilities/myvtkPolyDataNormals.cxx
#   matching/CreateDeformationField.cxx
# )

SET(VTK_LIBRARIES
  vtkChartsCore
  vtkRenderingLOD
  vtkRenderingOpenGL
  vtkRenderingVolume
  vtkRenderingVolumeOpenGL
  vtkViewsCore
  vtkViewsContext2D
  vtkInteractionStyle
  vtkIOExport
  vtkIOGeometry
  vtkIOLegacy
  vtkIOImage
  vtkIOCore
  vtkImagingCore
  vtkRenderingFreeTypeOpenGL
  vtkRenderingFreeType
  vtkRenderingCore
  vtkFiltersModeling
  vtkFiltersExtraction
  vtkFiltersSources
  vtkFiltersGeometry
  vtkFiltersGeneral
  vtkFiltersCore
  vtkCommonExecutionModel
  vtkCommonDataModel
  vtkCommonMisc
  vtkCommonSystem
  vtkCommonTransforms
  vtkCommonMath
  vtkCommonCore
  vtksys
  vtkftgl
  vtkfreetype
  vtkzlib
)

#
# Target link :
#
Set ( VTK_DEPENDANCIES "" )
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        vtkCommon
        vtkGraphics
        vtkIO
        vtkRendering
        vtkFiltering
    )
Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        ${VTK_LIBRARIES}
    )
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )

Target_Link_Libraries(sparseGeodesicRegression3CMRep
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
  ${PDE_LIBS}
  ${SYSTEM_LIBS}
  toms611
)

#Target_Link_Libraries(sparseGeodesicRegression3CMRep
#  ${ARMADILLO_LIBRARIES}
#  ${ITK_LIBRARIES}
#  ${VTK_DEPENDANCIES}
#  ${CUDA_LIBS}
#  ${OPENGL_LIBRARIES}
#  ${OPENGL_glu_LIBRARY}
#  ${SYSTEM_LIBS}
#  ${PDE_LIBS}
#)

#
# Dependencies (only for VTK < 6) :
#
If ( "${VTK_MAJOR_VERSION}" LESS 6 )

Add_Dependencies (sparseGeodesicRegression3CMRep
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )



# add tests
enable_testing()

add_test(NAME surface_matching WORKING_DIRECTORY ../../examples/surface_matching/ COMMAND sparseMatching3 paramDiffeos.xml paramSurface.xml sourceSurfaceMesh.vtk targetSurfaceMesh.vtk)





