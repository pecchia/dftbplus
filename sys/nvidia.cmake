#
# Toolchain file for
#
# NVIDIA/PGI Fortran compiler (tested with version 21.2) / GNU C-compiler
#
# Notes:
#
#  * CAVE PGI: The NVIDIA/PGI Fortran compiler has a long-standing tradition of being horribly
#    broken (as of 2021). Usually, there are always several (valid) Fortran constructs it can not
#    cope with, and also those, which you could get compiled yesterday, may be broken with their
#    actual compiler of today. Therefore, be prepared for false compiler errors, internal compiler
#    errors and segfaulting executables. But, there is of course always the chance, that this day
#    happens to be your PGI-lucky day. Challenge your luck! Start the build, cross your fingers,
#    and hope for the best...
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#


#
# Fortran compiler settings
#
set(GPU_FLAGS "-acc -gpu=cc80 -Mcuda -Mcudalib=cublas -lcusolver -Mallocatable=03 -Mbackslash" CACHE STRING "GPU Flags")
#set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GPU_FLAGS} -Mallocatable=03 -Mbackslash"
#  CACHE STRING "Build type independent Fortran compiler flags")
#set(GPU_FLAGS "-lcusolver -lcublas" CACHE STRING "GPU Flags")

set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GPU_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O2"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g"
  CACHE STRING "Fortran compiler flags for Debug build")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

set(FYPP_FLAGS "" CACHE STRING "Fypp preprocessor flags")


#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_DEBUG "-g"
  CACHE STRING "C compiler flags for Debug build")


#
# External libraries
#

# NOTE: Libraries with CMake export files (e.g. ELSI and if the HYBRID_CONFIG_METHODS variable
# contains the "Find" method also libNEGF, libMBD, ScalapackFx and MpiFx) are included by searching
# for the export file in the paths defined in the CMAKE_PREFIX_PATH **environment** variable. Make
# sure your CMAKE_PREFIX_PATH variable is set up accordingly.

# LAPACK and BLAS
#set(LAPACK_LIBRARY "openblas" CACHE STRING "LAPACK and BLAS libraries to link")
#set(LAPACK_LIBRARY_DIR "" CACHE STRING
#  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
#set(ARPACK_LIBRARY "arpack" CACHE STRING "Arpack libraries")
#set(ARPACK_LIBRARY_DIR "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
#set(SCALAPACK_LIBRARY "scalapack-openmpi" CACHE STRING "Scalapack libraries to link")
#set(SCALAPACK_LIBRARY_DIR "" CACHE STRING "Directories where Scalapack libraries can be found")

# NOTE: The libraries below provide Pkg-Conf export files.  If your PKG_CONFIG_PATH environment
# variable has been set up correctly (containing the paths to these libraries), no adjustment should
# be necessary below.

# PLUMED -- only needed when compiled with PLUMED support
#set(PLUMED_LIBRARY "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")
#set(PLUMED_LIBRARY_DIR "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
#set(MAGMA_LIBRARY "magma" CACHE STRING "Magma library")
#set(MAGMA_LIBRARY_DIR "" CACHE STRING "Directories to scan for MAGMA library")
#set(MAGMA_INCLUDE_DIRECTORY "" CACHE STRING "Directories to scan for MAGMA include files")
