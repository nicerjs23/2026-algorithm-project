cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp
)

# RB-tree + DP 알고리즘
add_executable(redblackTree_DP
        redblackTree_DP.cpp)
target_link_libraries(redblackTree_DP psapi)

# bruteforce 알고리즘
add_executable(bruteforce
        bruteforce.cpp)
target_link_libraries(bruteforce psapi)

