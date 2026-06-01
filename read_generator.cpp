// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할2 (김세훈)
//  Read 생성기
//
//  reference_builder 가 만든 원본 파일을 자동 감지해 두 세트의 read 를 생성한다.
//   A) original_synthetic_1M.txt -> reads_synthetic.txt   (항상 생성)
//   B) original_yeast_1M.txt     -> reads_yeast.txt       (파일 있을 때만)
//
//  [설계 방향]
//   - read 자체에는 에러/변이를 넣지 않는다.
//   - 돌연변이는 reference 쪽에만 존재 (역할1 코드에서 처리됨)
//   - 알고리즘 팀(3,4번)은 reads_synthetic.txt 를 기본 입력으로 사용
//     빵효모 비교 실험 시 reads_yeast.txt 를 추가로 사용
//
//  [출력 형식]
//   read_id \t start_pos \t sequence
//   - read_id  : 0-based 인덱스
//   - start_pos: 해당 original 기준 0-based 시작 위치 (정답)
//   - sequence : READ_LENGTH 길이의 순수 ATCG 서열
//
//  [팀원 참고]
//   - SEED=42 고정 -> 팀원 누구나 동일한 reads 파일 생성
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <random>

// ──────────────── 파라미터 ────────────────
static const std::string IN_SYNTH_ORIGINAL  = "original_synthetic_1M.txt";
static const std::string IN_YEAST_ORIGINAL  = "original_yeast_1M.txt";

static const std::string OUT_SYNTH_READS    = "reads_synthetic.txt";
static const std::string OUT_YEAST_READS    = "reads_yeast.txt";

static const int      READ_LENGTH = 30;      // read 길이 (bp) - 팀 공통 L=30
static const int      NUM_READS   = 100000;  // 생성할 read 수 - 팀 공통 M=100,000
static const unsigned SEED        = 42;      // 고정 시드 (팀 공통)
// ──────────────────────────────────────────

// 주어진 원본 파일에서 랜덤 위치 read 를 NUM_READS 개 뽑아 out_reads 에 저장한다.
// 파일이 없으면 false 반환 (호출자가 옵션/필수 여부 판단).
static bool generateOne(const std::string &in_original,
                        const std::string &out_reads,
                        const std::string &label) {
    std::ifstream fin(in_original);
    if (!fin.is_open()) {
        return false;  // 호출자가 처리
    }
    std::string original;
    std::getline(fin, original);
    fin.close();

    if ((int)original.size() < READ_LENGTH) {
        std::cerr << "[오류] " << in_original << " 길이(" << original.size()
                  << ")가 READ_LENGTH(" << READ_LENGTH << ")보다 짧습니다.\n";
        return false;
    }
    std::cout << "       [" << label << "] " << in_original
              << " 로드 (" << original.size() << " 염기)\n";

    std::ofstream fout(out_reads);
    if (!fout.is_open()) {
        std::cerr << "[오류] " << out_reads << " 생성 실패.\n";
        return false;
    }
    fout << "read_id\tstart_pos\tsequence\n";

    std::mt19937 rng(SEED);
    std::uniform_int_distribution<int> startDist(0, (int)original.size() - READ_LENGTH);

    for (int i = 0; i < NUM_READS; ++i) {
        int start = startDist(rng);
        fout << i << '\t' << start << '\t'
             << original.substr(start, READ_LENGTH) << '\n';
    }

    std::cout << "       [" << label << "] " << out_reads
              << " 생성 완료 (read " << NUM_READS
              << "개 / L=" << READ_LENGTH << ")\n";
    return true;
}

bool generate_reads() {
    // ========== A. 인공 서열 read (항상 생성) ==========
    std::cout << "[A] 인공 서열 read 생성\n";
    if (!generateOne(IN_SYNTH_ORIGINAL, OUT_SYNTH_READS, "synthetic")) {
        std::cerr << "[오류] 인공 서열 read 생성 실패. "
                  << IN_SYNTH_ORIGINAL << " 가 같은 폴더에 있는지 확인하세요.\n";
        return false;
    }

    // ========== B. 빵효모 read (원본 파일 있을 때만) ==========
    {
        std::ifstream test(IN_YEAST_ORIGINAL);
        if (test.is_open()) {
            test.close();
            std::cout << "\n[B] " << IN_YEAST_ORIGINAL
                      << " 발견 -> 빵효모 read 도 함께 생성\n";
            if (!generateOne(IN_YEAST_ORIGINAL, OUT_YEAST_READS, "yeast")) {
                return false;
            }
        } else {
            std::cout << "\n[알림] " << IN_YEAST_ORIGINAL
                      << " 없음 -> 빵효모 read 건너뜀.\n";
        }
    }

    std::cout << "\n[완료] read 생성 단계 종료.\n"
              << "       알고리즘 팀은 reads_synthetic.txt 를 기본 입력으로,\n"
              << "       비교 실험 시 reads_yeast.txt 도 함께 사용하세요.\n";
    return true;
}
