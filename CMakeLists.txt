cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

add_executable(2026_algorithm_project main.cpp)

# 역할2 (김세훈): Read 생성기
add_executable(read_generator read_generator.cpp)

# 역할2 (김세훈): BWT + FM-index 매퍼
add_executable(bwt_mapper bwt_mapper.cpp)
