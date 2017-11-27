cmake_minimum_required(VERSION 3.8)
project(NanotubeBuilder)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set(SOURCE_FILES src/main.cpp
        src/NanotubeBuilder.h src/NanotubeBuilder.cpp
        src/LatticeOperations.h src/LatticeOperations.cpp src/PMath.h)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output)

add_subdirectory(tests)

add_executable(NanotubeBuilder ${SOURCE_FILES})