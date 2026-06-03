cmake_minimum_required(VERSION 4.0)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# 파이프라인: main -> reference_builder -> read_generator
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)

# RB-tree + DP 알고리즘
add_executable(redblackTree_DP
        redblackTree_DP.cpp)

# bruteforce 알고리즘
add_executable(bruteforce
        bruteforce.cpp)

# Hash Table + Greedy 알고리즘 (역할3, 이동건)
add_executable(hash_greedy
        hash_greedy.cpp)

# Hash Table + DP 알고리즘 (역할3, 이동건) - 인델 환경 비교용
add_executable(hash_dp
        hash_dp.cpp)

# Windows 에서만 psapi 라이브러리 링크 (메모리 측정용)
# macOS/Linux 는 POSIX getrusage 를 쓰므로 추가 라이브러리 불필요
if(WIN32)
    target_link_libraries(redblackTree_DP psapi)
    target_link_libraries(bruteforce      psapi)
    target_link_libraries(hash_greedy     psapi)
    target_link_libraries(hash_dp         psapi)
endif()