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

AC_INIT(Bio++ Utils, 1.3.1, Julien.Dutheil@univ-montp2.fr, bpp-utils)
AC_CONFIG_SRCDIR(Utils/Clonable.h)
AM_INIT_AUTOMAKE()
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_MSG_CHECKING([Covariant return with virtual inheritance])
AC_ARG_WITH(
	[virtual-cov],
	[  --with-virtual-cov      Enable covariant return type with
	                  virtual inheritance (GCC >= 3.4)],
  ,
  with_virtual_cov="yes"
)
AC_MSG_RESULT($with_virtual_cov);
if test "$with_virtual_cov" == "no"
then
	AC_MSG_WARN([
+--------------------------------------------------------+
|Covariant return type with virtual inheritance disabled.|
+--------------------------------------------------------+
]);
	AC_DEFINE(NO_VIRTUAL_COV, yes)
fi
AC_OUTPUT(Makefile Utils/Makefile)
