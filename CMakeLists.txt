cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)

# 역할2 (김세훈): Trap Read 생성기 (Baseline / InDel / End-Heavy)
add_executable(trap_generator trap_generator.cpp)

# 역할2 (김세훈): BWT + FM-index 매퍼
add_executable(bwt_mapper bwt_mapper.cpp)
