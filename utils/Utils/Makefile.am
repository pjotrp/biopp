#Copyright or © or Copr. CNRS, (November 16, 2004-2008)
#
#This software is a computer program whose purpose is to provide utilitary
#classes. This file belongs to the Bio++ Project.
#
#This software is governed by the CeCILL  license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL license and that you accept its terms.

H_FILES = \
StringTokenizer.h\
Exceptions.h\
MapTools.h\
AttributesTools.h\
TextTools.h\
FileTools.h\
Number.h\
Clonable.h\
IOFormat.h\
ApplicationTools.h\
exceptions\
BppString.h\
BppVector.h\
RGBColor.h\
ColorTools.h\
ColorManager.h\
ColorSet.h\
DefaultColorSet.h\
RColorSet.h\
MolscriptColorSet.h\
DvipsColorSet.h\
Font.h\
FontManager.h\
GraphicDevice.h\
XFigGraphicDevice.h\
SVGGraphicDevice.h\
PGFGraphicDevice.h\
KeyvalTools.h

CPP_FILES = \
Exceptions.cpp\
AttributesTools.cpp\
TextTools.cpp\
FileTools.cpp\
StringTokenizer.cpp\
ApplicationTools.cpp\
ColorTools.cpp\
RColorSet.cpp\
MolscriptColorSet.cpp\
DvipsColorSet.cpp\
GraphicDevice.cpp\
XFigGraphicDevice.cpp\
SVGGraphicDevice.cpp\
PGFGraphicDevice.cpp\
KeyvalTools.cpp


lib_LTLIBRARIES = libbpp-utils.la
libbpp_utils_la_SOURCES = $(H_FILES) $(CPP_FILES)
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current: The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age: The difference between the newest and oldest interfaces that this library implements.
# In other words, the library implements all the interface numbers in the range from number current - age to current.
libbpp_utils_la_LDFLAGS = -version-info 5:1:1

library_includedir=$(includedir)/Utils
library_include_HEADERS = $(H_FILES)
