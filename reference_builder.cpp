// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할1 (이동건)
//  Reference 게놈 생성기
//
//  main.cpp 에서 SNP 비율 하나만 넘기면
//   1) FASTA -> original_1M.txt        (원본 100만 염기)
//   2) original -> reference_genome.txt (SNP 삽입본)
//                + snp_list.txt          (변이 위치 정답표)
//  까지 한 번에 처리한다.
//
//  실행 위치: CLion 은 기본적으로 cmake-build-debug/ 에서 실행되며,
//             sequence.fasta 도 그 폴더에 들어있으므로 결과 파일도 같은 폴더에 생성된다.
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <array>

// ---------------- 파라미터 (필요시 여기만 수정) ----------------
static const std::string INPUT_FASTA   = "sequence.fasta";        // 다운받은 원본 FASTA
static const std::string OUT_ORIGINAL  = "original_1M.txt";       // 1단계 결과: 원본 100만
static const std::string OUT_REFERENCE = "reference_genome.txt";  // 2단계 결과: SNP 삽입본
static const std::string OUT_SNP_LIST  = "snp_list.txt";          // 변이 위치 정답표

static const std::size_t N    = 1000000;  // 원본 서열 길이 (100만)
static const unsigned    SEED = 42;       // 고정 시드 -> 팀원 모두 동일 데이터 생성
// -------------------------------------------------------------

// FASTA 에서 헤더(>로 시작)를 제외하고 ATCG 염기만 모아 N개 추출한다.
static std::string extractOriginal(const std::string &fastaPath, std::size_t need) {
    std::ifstream in(fastaPath);
    if (!in.is_open()) {
        std::cerr << "[오류] " << fastaPath << " 열기 실패. 같은 폴더에 있는지 확인하세요.\n";
        return "";
    }

    std::string seq;
    seq.reserve(need);

    std::string line;
    while (std::getline(in, line) && seq.size() < need) {
        if (!line.empty() && line[0] == '>') continue;  // FASTA 헤더 줄 건너뛰기
        for (char c : line) {
            if (c == 'a' || c == 'A') c = 'A';
            else if (c == 'c' || c == 'C') c = 'C';
            else if (c == 'g' || c == 'G') c = 'G';
            else if (c == 't' || c == 'T') c = 'T';
            else continue;  // 공백, \r, N 등 ATCG 가 아닌 문자는 무시
            seq.push_back(c);
            if (seq.size() >= need) break;
        }
    }
    return seq;
}

// 염기 한 글자를 자기 자신이 아닌 다른 염기로 바꾼다 (= SNP).
static char mutateBase(char base, std::mt19937 &rng) {
    static const std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};
    int idx;
    switch (base) {
        case 'A': idx = 0; break;
        case 'C': idx = 1; break;
        case 'G': idx = 2; break;
        case 'T': idx = 3; break;
        default:  return base;  // ATCG 가 아니면 변이 안 함 (안전장치)
    }
    std::uniform_int_distribution<int> offset(1, 3);  // 1~3 더해서 mod 4 -> 항상 다른 염기
    return BASES[(idx + offset(rng)) & 3];
}

// SNP 비율(예: 0.001 = 0.1%)을 받아 원본/레퍼런스/SNP 목록 파일을 생성한다.
bool build_reference(double snp_rate) {
    // ---------- 1단계: FASTA -> 원본 100만 염기 ----------
    std::cout << "[1단계] " << INPUT_FASTA << " 에서 원본 " << N << " 염기 추출 중...\n";
    std::string original = extractOriginal(INPUT_FASTA, N);

    if (original.size() < N) {
        std::cerr << "[오류] 추출된 염기 수가 부족합니다. (" << original.size()
                  << " / " << N << ") FASTA 파일이 충분히 긴지 확인하세요.\n";
        return false;
    }

    {
        std::ofstream out(OUT_ORIGINAL);
        if (!out.is_open()) {
            std::cerr << "[오류] " << OUT_ORIGINAL << " 생성 실패.\n";
            return false;
        }
        out << original;  // 헤더 없이 서열만 한 줄로 저장
    }
    std::cout << "       -> " << OUT_ORIGINAL << " 생성 완료 (" << original.size() << " 염기)\n";

    // ---------- 2단계: 원본 -> SNP 삽입 -> reference 게놈 ----------
    std::cout << "[2단계] SNP " << (snp_rate * 100) << "% 삽입하여 reference 게놈 생성 중...\n";
    std::string reference = original;  // 원본을 복사한 뒤 일부 글자만 변이시킨다

    std::mt19937 rng(SEED);
    std::uniform_real_distribution<double> prob(0.0, 1.0);

    std::ofstream snpOut(OUT_SNP_LIST);
    if (!snpOut.is_open()) {
        std::cerr << "[오류] " << OUT_SNP_LIST << " 생성 실패.\n";
        return false;
    }
    snpOut << "position\toriginal\treference\n";  // 0-based 위치, 원본염기, 바뀐염기

    std::size_t snpCount = 0;
    for (std::size_t i = 0; i < reference.size(); ++i) {
        if (prob(rng) < snp_rate) {
            char before = reference[i];
            char after  = mutateBase(before, rng);
            reference[i] = after;
            snpOut << i << '\t' << before << '\t' << after << '\n';
            ++snpCount;
        }
    }

    {
        std::ofstream out(OUT_REFERENCE);
        if (!out.is_open()) {
            std::cerr << "[오류] " << OUT_REFERENCE << " 생성 실패.\n";
            return false;
        }
        out << reference;
    }

    std::cout << "       -> " << OUT_REFERENCE << " 생성 완료\n";
    std::cout << "       -> " << OUT_SNP_LIST << " 생성 완료 (변이 위치 정답표)\n";

    // ---------- 결과 요약 ----------
    double actualRate = static_cast<double>(snpCount) / reference.size() * 100.0;
    std::cout << "[완료] 총 길이 " << reference.size()
              << " 중 SNP " << snpCount << "개 삽입 "
              << "(실제 비율 " << actualRate << "%)\n";
    std::cout << "      reference 게놈은 original 과 정확히 " << snpCount
              << " 글자만 다릅니다.\n\n";

    return true;
}
