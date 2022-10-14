cmake_minimum_required(VERSION 3.23)
project(SKIPOD_2)
include_directories(/opt/homebrew/Cellar/open-mpi/4.1.4_2/include)

find_package(MPI REQUIRED)
set(CMAKE_CXX_STANDARD 11)

add_executable(SKIPOD_2 main.cpp)
target_link_libraries(SKIPOD_2 ${MPI_C_LIBRARIES})

