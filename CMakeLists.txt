cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)
