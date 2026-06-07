cmake_minimum_required(VERSION 3.20)
project(2026_algorithm_project)

set(CMAKE_CXX_STANDARD 17)

# ===== 데이터 생성 파이프라인 (역할1·2) =====
# main 실행 시 reference + read 데이터 한 번에 생성
add_executable(2026_algorithm_project
        main.cpp
        reference_builder.cpp
        read_generator.cpp)

# 12개 trap read 데이터셋 생성기 (역할2)
add_executable(trap_generator
        trap_generator.cpp)

# SNP 비율별 데이터셋 일괄 생성기 (역할3, 이동건) - hash 매트릭스 실험용
add_executable(dataset_generator
        dataset_generator.cpp)

# ===== 매칭 알고리즘들 =====
# bruteforce (모든 알고리즘 비교의 벤치마크)
add_executable(bruteforce
        bruteforce.cpp)

# Hash Table + Greedy (역할3, 이동건)
add_executable(hash_greedy
        hash_greedy.cpp)

# Hash Table + DP (역할3, 이동건) - 인델 환경 비교용
add_executable(hash_dp
        hash_dp.cpp)

# RB-tree + DP (역할2)
add_executable(redblackTree_DP
        redblackTree_DP.cpp)

# Suffix Array + Minimizer (역할3)
add_executable(suffix_array
        Suffix_Array_Minimizer.cpp)

# BWT + Suffix Array (역할4)
add_executable(bwt_mapper
        bwt_mapper.cpp)

# ===== Windows 전용: psapi 링크 (메모리 측정용) =====
# macOS/Linux 는 POSIX getrusage 를 쓰므로 추가 라이브러리 불필요
if(WIN32)
    target_link_libraries(bruteforce      psapi)
    target_link_libraries(hash_greedy     psapi)
    target_link_libraries(hash_dp         psapi)
    target_link_libraries(redblackTree_DP psapi)
endif()
