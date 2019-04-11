#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.11)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
#use intel compiler
set (CMAKE_CXX_COMPILER "mpiicpc")
set (CMAKE_C_COMPILER "mpiicc")

#set as inlcude dir the current directory
#one could devide the project into more directorys later
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(project_sources
	MY_MCP.cpp
	costs.cpp
	utils.cpp
        )

#add compiler flags inlcuding cblas
set(C_flags	"-mkl -Ofast -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS	"${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project (MY_MCP)
#add an executable - using variables
add_executable(${PROJECT_NAME} ${project_sources})

#include MKL - here include your mkl path
include_directories("/opt/intel/mkl/include")
#place where libfabric is installed - needed for IntelMPI-2019
target_link_libraries(${PROJECT_NAME} "/usr/local/lib/libfabric.so.1")
