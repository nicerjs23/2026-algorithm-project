// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할1 (이동건)
//  Reference 게놈 생성기
//
//  main.cpp 에서 SNP 비율 하나만 넘기면 두 종류의 원본을 만든다.
//   A) 인공 원본 (mt19937, ATCG 25% 균등 랜덤)         -- 항상 생성 (메인 baseline)
//      -> original_synthetic_1M.txt
//      -> reference_synthetic.txt + snp_list_synthetic.txt
//   B) 빵효모 원본 (sequence.fasta 에서 추출)            -- 파일 있을 때만 생성 (옵션)
//      -> original_yeast_1M.txt
//      -> reference_yeast.txt + snp_list_yeast.txt
//
//  read_generator 가 두 원본을 자동 감지해 reads_synthetic.txt / reads_yeast.txt
//  로 각각 read 를 뽑으므로 별도 별칭 파일은 필요 없다.
//
//  실행 위치: CLion 은 cmake-build-debug/ 에서 실행되며 결과도 같은 폴더에 생성된다.
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <array>

// ---------------- 파라미터 ----------------
static const std::string INPUT_FASTA = "sequence.fasta";  // 빵효모 (옵션)

// 인공 서열 결과 파일 (메인)
static const std::string OUT_SYNTH_ORIGINAL  = "original_synthetic_1M.txt";
static const std::string OUT_SYNTH_REFERENCE = "reference_synthetic.txt";
static const std::string OUT_SYNTH_SNP_LIST  = "snp_list_synthetic.txt";

// 빵효모 결과 파일 (sequence.fasta 있을 때만)
static const std::string OUT_YEAST_ORIGINAL  = "original_yeast_1M.txt";
static const std::string OUT_YEAST_REFERENCE = "reference_yeast.txt";
static const std::string OUT_YEAST_SNP_LIST  = "snp_list_yeast.txt";

static const std::size_t N    = 1000000;  // 원본 서열 길이 (100만)
static const unsigned    SEED = 42;       // 고정 시드 -> 팀원 모두 동일 데이터 생성
// -------------------------------------------------------------

// FASTA 에서 헤더(>로 시작)를 제외하고 ATCG 염기만 모아 need 개 추출한다.
static std::string extractOriginal(const std::string &fastaPath, std::size_t need) {
    std::ifstream in(fastaPath);
    if (!in.is_open()) return "";

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

// 난수 기반 인공 원본 서열 생성 (ATCG 25% 균등).
// 빵효모 같은 실제 서열은 반복(repeat)/편향 GC content 때문에
// 알고리즘의 순수 성능을 비교하기에 적합하지 않다.
// -> mt19937 로 균등 랜덤 ATCG 를 뽑아 통제된 baseline 을 만든다.
static std::string generateSynthetic(std::size_t need, std::mt19937 &rng) {
    static const std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> pick(0, 3);

    std::string seq;
    seq.reserve(need);
    for (std::size_t i = 0; i < need; ++i) {
        seq.push_back(BASES[pick(rng)]);
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

// 주어진 원본에 SNP 를 심고 reference / snp_list 파일을 저장한다.
// label 은 로그 출력용 ("synthetic" / "yeast").
static bool mutateAndSave(const std::string &original,
                          double snp_rate,
                          unsigned seed,
                          const std::string &out_reference,
                          const std::string &out_snp_list,
                          const std::string &label) {
    std::string reference = original;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> prob(0.0, 1.0);

    std::ofstream snpOut(out_snp_list);
    if (!snpOut.is_open()) {
        std::cerr << "[오류] " << out_snp_list << " 생성 실패.\n";
        return false;
    }
    snpOut << "position\toriginal\treference\n";

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

    std::ofstream out(out_reference);
    if (!out.is_open()) {
        std::cerr << "[오류] " << out_reference << " 생성 실패.\n";
        return false;
    }
    out << reference;

    double actualRate = static_cast<double>(snpCount) / reference.size() * 100.0;
    std::cout << "       [" << label << "] " << out_reference
              << " 생성 (SNP " << snpCount << "개, 실제 " << actualRate << "%)\n";
    std::cout << "       [" << label << "] " << out_snp_list << " 생성 (정답표)\n";
    return true;
}

// SNP 비율(예: 0.001 = 0.1%)을 받아 인공/빵효모 원본+레퍼런스+SNP 목록을 만든다.
bool build_reference(double snp_rate) {
    // ========== A. 인공 서열 (메인 baseline, 항상 생성) ==========
    std::cout << "[A] 인공 원본 서열 생성 (mt19937, ATCG 25% 균등, " << N << " 염기)\n";
    std::mt19937 synth_rng(SEED);
    std::string synthetic = generateSynthetic(N, synth_rng);

    {
        std::ofstream out(OUT_SYNTH_ORIGINAL);
        if (!out.is_open()) {
            std::cerr << "[오류] " << OUT_SYNTH_ORIGINAL << " 생성 실패.\n";
            return false;
        }
        out << synthetic;
    }
    std::cout << "       [synthetic] " << OUT_SYNTH_ORIGINAL << " 생성 완료\n";

    std::cout << "       SNP " << (snp_rate * 100) << "% 삽입 중...\n";
    if (!mutateAndSave(synthetic, snp_rate, SEED,
                       OUT_SYNTH_REFERENCE, OUT_SYNTH_SNP_LIST, "synthetic")) {
        return false;
    }

    // ========== B. 빵효모 (sequence.fasta 있을 때만) ==========
    std::ifstream test(INPUT_FASTA);
    if (test.is_open()) {
        test.close();
        std::cout << "\n[B] " << INPUT_FASTA << " 발견 -> 빵효모 서열도 함께 처리\n";

        std::string yeast = extractOriginal(INPUT_FASTA, N);
        if (yeast.size() < N) {
            std::cerr << "[경고] 빵효모 추출 부족 (" << yeast.size() << "/" << N
                      << ") -> 빵효모 처리 건너뜀.\n";
        } else {
            std::ofstream out(OUT_YEAST_ORIGINAL);
            if (!out.is_open()) {
                std::cerr << "[오류] " << OUT_YEAST_ORIGINAL << " 생성 실패.\n";
                return false;
            }
            out << yeast;
            std::cout << "       [yeast] " << OUT_YEAST_ORIGINAL
                      << " 생성 완료 (" << yeast.size() << " 염기)\n";

            std::cout << "       SNP " << (snp_rate * 100) << "% 삽입 중...\n";
            if (!mutateAndSave(yeast, snp_rate, SEED,
                               OUT_YEAST_REFERENCE, OUT_YEAST_SNP_LIST, "yeast")) {
                return false;
            }
        }
    } else {
        std::cout << "\n[알림] " << INPUT_FASTA
                  << " 없음 -> 빵효모 처리 건너뜀 (인공 서열만 사용).\n";
    }

    std::cout << "\n[완료] reference 생성 단계 종료.\n\n";
    return true;
}
