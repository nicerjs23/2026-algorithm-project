cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)

# 역할3 (권혁윤): 브루트포스 벤치마크
add_executable(brute_force brute_force.cpp)
target_link_libraries(brute_force psapi)