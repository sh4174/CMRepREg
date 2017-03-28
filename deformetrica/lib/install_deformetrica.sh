#!/bin/bash

# =========================================================================== #
#                    Automatic installation of Deformetrica                   #
# Description :                                                               #
# -------------                                                               #
#   This script will automatically install the CMake software together with   #
# the ITK & VTK libraries and compile the project. For the moment,            #
#                                                                             #
# TODO :                                                                      #
# ------                                                                      #
# * check if CMake is up-to-date (current version : automatic                 #
# installation of CMake)                                                      #
# =========================================================================== #


# Path to download the several softwares :
path_to_download_cmake="http://www.cmake.org/files/v2.8/cmake-2.8.10.2.tar.gz"
path_to_download_itk="http://sourceforge.net/projects/itk/files/itk/4.3/InsightToolkit-4.3.1.tar.gz/download"
path_to_download_vtk="http://www.vtk.org/files/release/5.10/vtk-5.10.1.tar.gz"


# Note : " > /dev/null " is used to quiet messages when executing commands

path_to_make="/usr/bin/make"


echo "##################################################"
echo "# Installation of cmake :"
echo "# (if not installed or not up-to-date)"
echo "##################################################"
path_to_cmake="/lena16/alzheimer1/Software/CMake/cmakeb--2.8.10.1/bin/cmake"
#path_to_cmake=`which cmake`
#if [ $? != 0]; then # No 'cmake' command found
	# Manual installation of CMake :
	# TODO
	wget "$path_to_download_cmake" -O cmake-latest.tar.gz
	echo "Extracting CMake..."
	tar -zxf cmake-latest.tar.gz
	rm -v cmake-latest.tar.gz
	folder_source_cmake=`ls | grep cmake-[0-9]*`
	if [ -d $folder_source_cmake ]; then
		echo "Installation of $folder_source_cmake"
		mkdir cmake
		path_to_cmake=`pwd`"/cmake/bin/cmake"
		cd $folder_source_cmake
		./bootstrap --prefix=../cmake && $path_to_make && $path_to_make install
	else
		echo "Error while extracting CMake"
		exit 1
	fi

#	echo "CMake program not found or not up-to-date."
#	exit 1
#fi

cd ../
rm -rf $folder_source_cmake

echo "##################################################"
echo "# Installation of ITK :"
echo "##################################################"
read -p 'Enter the path to ITK (if not installed, press ENTER)  : ' path_to_itk
if [ -z $path_to_itk ]; then
	wget "$path_to_download_itk" -O InsightToolkit-latest.tar.gz
	if [ $? != 0 ]; then
		echo "Error when downloading ITK"
		exit 1
	else
		# Extracting ITK :
		echo "Extracting ITK..."
		tar -zxf InsightToolkit-latest.tar.gz
		rm -v InsightToolkit-latest.tar.gz
		folder_source_itk=`ls | grep InsightToolkit-[0-9]*`
		if [ -d $folder_source_itk ]; then
			# Installation of ITK :
			mkdir ITKb;
			cd ITKb/;
			$path_to_cmake -D USE_FFTWD=ON -D USE_FFTWF=ON ../$folder_source_itk/
			$path_to_make
			$path_to_itk=`pwd`
			# Cleaning :
			cd ../
			rm -rf $folder_source_itk
		else
			echo "Error while extracting ITK"
			exit 1
		fi
	fi
fi


echo "##################################################"
echo "# Installation of VTK :"
echo "##################################################"
read -p 'Enter the path to VTK (if not installed, press ENTER)  : ' path_to_vtk
if [ -z $path_to_vtk ]; then
	wget "$path_to_download_vtk" -O vtk-latest.tar.gz
	if [ $? != 0 ]; then
		echo "Error when downloading VTK."
		exit 1
	else
		echo "Extracting VTK..."
		tar -zxf vtk-latest.tar.gz
		rm -v vtk-latest.tar.gz
		folder_source_vtk=`ls | grep VTK[0-9]*`
		if [ -d $folder_source_vtk ]; then
			# Installation of VTK :
			mkdir VTKb
			cd VTKb/
			$path_to_cmake ../$folder_source_vtk/
			$path_to_make
			$path_to_vtk=`pwd`
			# Cleaning :
			cd ..
			rm -rf $folder_source_vtk
		else
			echo "Error while extracting the VTK."
			exit 1
		fi
	fi
fi

mkdir ../bin 2> /dev/null 
cd ../bin/


echo "##################################################"
echo "# Compilation of deformatrica :"
echo "##################################################"
$path_to_cmake -D CMAKE_BUILD_TYPE=Release -D ITK_DIR=$path_to_itk -D VTK_DIR=$path_to_vtk ../app
$path_to_make

