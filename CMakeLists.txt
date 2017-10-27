cmake_minimum_required(VERSION 3.8)
project(NanotubeBuilder)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/main.cpp src/NanotubeBuilder.h src/NanotubeBuilder.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output)
add_subdirectory(tests)
include_directories(src)
add_executable(NanotubeBuilder ${SOURCE_FILES})