// =============================================================
//  Short-read re-sequencing 프로젝트 - 파이프라인 진입점
//
//  여기서는 단계별 함수만 순서대로 호출한다.
//  실제 로직은 각 .cpp 에 분리되어 있다.
//
// =============================================================

#include <iostream>

// Windows 콘솔에서 한글 깨짐 방지용. macOS/Linux 는 UTF-8 이 기본이라 불필요.
#ifdef _WIN32
#include <windows.h>
#endif

// 다른 .cpp 파일에 정의된 함수 선언 (헤더 대신 한 줄로)
bool build_reference(double snp_rate);
bool generate_reads();

// ---------------- 파라미터 ----------------
static const double SNP_RATE = 0.001;  // SNP 비율: 0.001 = 0.1%
// -----------------------------------------

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8 (한글 깨짐 방지)
#endif

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
