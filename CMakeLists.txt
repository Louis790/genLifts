cmake_minimum_required(VERSION 3.27)
project(genLifts)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src/graph)
add_subdirectory(src/group)
add_subdirectory(src/algorithms)
add_subdirectory(src/nautyAndMultigraph)

add_executable(genLifts
         src/Main.cpp
)

add_executable(dreadnaut
         src/nautyAndMultigraph/dreadnaut.c
)

target_link_libraries(genLifts PUBLIC Graph Group Algorithms nautyAndMultigraph)
target_link_libraries(dreadnaut PUBLIC nautyAndMultigraph)
