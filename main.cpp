// =============================================================
//  Short-read re-sequencing 프로젝트 - 파이프라인 진입점
//
//  여기서는 단계별 함수만 순서대로 호출한다.
//  실제 로직은 각 .cpp 에 분리되어 있다.
//
//   1단계 (역할1, 이동건): reference_builder.cpp
//        FASTA -> original_1M.txt
//                + reference_genome.txt (SNP 삽입본)
//                + snp_list.txt          (변이 위치 정답표)
//
//   2단계 (역할2, 김세훈): read_generator.cpp
//        original_1M.txt -> reads.txt
//
//  SNP 비율은 여기 main 에서 SNP_RATE 상수 하나만 바꾸면 된다.
//   - 0.001 (0.1%) / 0.005 (0.5%) / 0.01 (1%) 세 가지를 돌려가며 실험.
// =============================================================

#include <iostream>

// 다른 .cpp 파일에 정의된 함수 선언 (헤더 대신 한 줄로)
bool build_reference(double snp_rate);
bool generate_reads();

// ---------------- 파라미터 ----------------
static const double SNP_RATE = 0.05;  // SNP 비율: 0.001 = 0.1%
// -----------------------------------------

int main() {
    std::cout << "=== [1단계] reference 게놈 생성 (SNP "
              << (SNP_RATE * 100) << "%) ===\n";
    if (!build_reference(SNP_RATE)) {
        std::cerr << "[중단] reference 생성 실패.\n";
        return 1;
    }

    std::cout << "=== [2단계] read 생성 ===\n";
    if (!generate_reads()) {
        std::cerr << "[중단] read 생성 실패.\n";
        return 1;
    }

    std::cout << "\n=== 파이프라인 완료 ===\n";
    return 0;
}