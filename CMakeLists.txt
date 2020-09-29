cmake_minimum_required(VERSION 3.8)
project(NanoBuilder)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
include_directories(src)
file(GLOB SOURCE_FILES "src/*.cpp")
file(GLOB HEADERS "src/*.h" )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/output)

add_subdirectory(tests)

add_executable(NanoBuilder ${SOURCE_FILES} ${HEADERS})