// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할2 (김세훈)
//  Read 생성기: original 게놈에서 read 깨끗하게 잘라내기
//
//  [설계 방향]
//   - read 자체에는 에러/변이를 넣지 않는다.
//   - 돌연변이 비율(0.1% / 0.5% / 1%)은 1번(이동건) 코드에서
//     reference_genome 을 3종 생성하는 방식으로 반영한다.
//   - 알고리즘(3번 KMP, 4번 Seed-and-Extend)은 깨끗한 read 를
//     각 비율의 reference_genome 에서 탐색한다.
//
//  [파이프라인]
//   original_1M.txt → (랜덤 위치에서 100bp 잘라내기) → reads.txt
//
//  [출력 형식 - reads.txt]
//   read_id \t start_pos \t sequence
//   - read_id  : 0-based 인덱스
//   - start_pos: original_1M.txt 기준 0-based 시작 위치 (정답)
//   - sequence : READ_LENGTH 길이의 순수 ATCG 서열
//
//  [팀원 참고]
//   - 3번(KMP), 4번(Seed-and-Extend) 모두 reads.txt 를 공통 입력으로 사용
//   - SEED=42 고정 → 팀원 누구나 동일한 reads.txt 생성 가능
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

// ──────────────── 파라미터 (필요시 여기만 수정) ────────────────
static const std::string IN_ORIGINAL  = "original_1M.txt";
static const std::string OUT_READS    = "reads.txt";

static const int      READ_LENGTH = 30;      // read 길이 (bp) - 팀 공통 L=30
static const int      NUM_READS   = 100000;  // 생성할 read 수 - 팀 공통 M=100,000
static const unsigned SEED        = 42;     // 고정 시드 (팀 공통)
// ──────────────────────────────────────────────────────────────

// original_1M.txt 에서 랜덤 위치를 잡아 길이 READ_LENGTH 의 read 를 NUM_READS 개 생성한다.
// 성공하면 true, 실패하면 false 반환.
bool generate_reads() {
    // ── original 로드 ──
    std::ifstream fin(IN_ORIGINAL);
    if (!fin.is_open()) {
        std::cerr << "[오류] " << IN_ORIGINAL << " 열기 실패. 같은 폴더에 있는지 확인하세요.\n";
        return false;
    }
    std::string original;
    std::getline(fin, original);
    fin.close();

    if ((int)original.size() < READ_LENGTH) {
        std::cerr << "[오류] original 길이(" << original.size()
                  << ")가 READ_LENGTH(" << READ_LENGTH << ")보다 짧습니다.\n";
        return false;
    }
    std::cout << "[로드] original 길이: " << original.size() << " 염기\n";

    // ── read 생성 ──
    std::ofstream fout(OUT_READS);
    if (!fout.is_open()) {
        std::cerr << "[오류] " << OUT_READS << " 생성 실패.\n";
        return false;
    }
    fout << "read_id\tstart_pos\tsequence\n";

    std::mt19937 rng(SEED);
    std::uniform_int_distribution<int> startDist(0, (int)original.size() - READ_LENGTH);

    for (int i = 0; i < NUM_READS; ++i) {
        int start = startDist(rng);
        std::string seq = original.substr(start, READ_LENGTH);
        fout << i << '\t' << start << '\t' << seq << '\n';
    }

    std::cout << "[완료] " << OUT_READS << " 생성 완료\n";
    std::cout << "       read " << NUM_READS << "개 / 길이 " << READ_LENGTH << "bp / SEED=" << SEED << "\n";
    std::cout << "       3번(KMP), 4번(Seed-and-Extend)은 이 파일을 공통 입력으로 사용하세요.\n";

    return true;
}
