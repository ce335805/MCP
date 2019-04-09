#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.11)

#set as inlcude dir the current directory
#one could devide the project into more directorys later
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(project_sources
	MY_MCP.cpp
	costs.cpp
	utils.cpp
        )

#add compiler flags inlcuding cblas
set(C_flags	"-lcblas -O3 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS	"${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project (MY_MCP)
#add an executable - using variables
add_executable(${PROJECT_NAME} ${project_sources})
