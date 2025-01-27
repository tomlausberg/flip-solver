cmake_minimum_required( VERSION 2.6 )


SET( FLAGGY_FLAGS "-g -Wall -O2 -fopenmp")
SET (CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${FLAGGY_FLAGS}")
project( flip_solver )


find_package (Eigen3 3.3 REQUIRED NO_MODULE)

file(GLOB bev_files
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp
)
add_executable(flip_solver ${bev_files})

target_link_libraries( flip_solver Eigen3::Eigen)
