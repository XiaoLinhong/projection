cmake_minimum_required(VERSION 3.5)

project(projection LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -Ofast")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")

get_directory_property(parent_dir PARENT_DIRECTORY)

set(sources
  src/module_proj.f90
  src/projection.F90
)

add_library(projection ${sources})

if (NOT parent_dir)
  add_executable(projection.exe app/main.f90)
  target_link_libraries(projection.exe projection)
endif ()
