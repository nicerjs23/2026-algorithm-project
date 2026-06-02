cmake_minimum_required(VERSION 4.2)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 20)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)

add_executable(trap_generator trap_generator.cpp)

# suffix array + minimizer 연동
add_executable(suffix_array Suffix_Array_Minimizer.cpp)
