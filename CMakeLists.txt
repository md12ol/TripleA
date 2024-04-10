cmake_minimum_required(VERSION 3.16.5)
project(TripleA)

set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_BUILD_TYPE "Release")

add_executable(TripleA
        main.cpp
        main.h
        SDA.h
        SDA.cpp)