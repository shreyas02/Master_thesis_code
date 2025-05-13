cmake_minimum_required(VERSION 3.20)
project(mpicpp NONE)

# Put runtime binaries (executables) in the shared folder at the project root
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")
# Put shared libraries in the same shared folder
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")
# Put static libraries (archives) in the same shared folder
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")

# Create the shared library target instead of an executable target.
add_library(mpicpp SHARED mpicpp.cpp)

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# Get Trilinos as one entity but require the packages being used
find_package(Trilinos REQUIRED COMPONENTS Amesos2 Teuchos Tpetra Belos ShyLU_DDFROSch Xpetra)

# Find MPI package
#find_package(MPI REQUIRED)

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

# Make sure to use same compilers and flags as Trilinos
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )
set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

enable_language(C)
enable_language(CXX)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()


# Set the C++ standard and require it
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

execute_process(
  COMMAND julia --project=/path/from/where/you/execute/ -e "using libcxxwrap_julia_jll; print(libcxxwrap_julia_jll.artifact_dir)"
  OUTPUT_VARIABLE LIBCXXWRAP_ARTIFACT_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Linking all directories for CxxWrap, Julia Headers and Trilinos 
target_include_directories(mpicpp PRIVATE
${LIBCXXWRAP_ARTIFACT_DIR}/include
/path/to/julia/iso
${MPI_CXX_INCLUDE_PATH}
${CMAKE_CURRENT_SOURCE_DIR}
${Trilinos_INCLUDE_DIRS}
${Trilinos_TPL_INCLUDE_DIRS}
)

# Link with the necessary libraries including trilinos 
target_link_libraries(mpicpp PRIVATE 
${LIBCXXWRAP_ARTIFACT_DIR}/include
/path/to/julia/iso
#MPI::MPI_CXX
${Trilinos_LIBRARIES} 
${Trilinos_TPL_LIBRARIES}
)

# MPI Related
target_compile_definitions(mpicpp PRIVATE ${MPI_CXX_COMPILE_DEFINITIONS})
target_compile_options(mpicpp PRIVATE ${MPI_CXX_COMPILE_OPTIONS})
target_link_options(mpicpp PRIVATE ${MPI_CXX_LINK_FLAGS})
