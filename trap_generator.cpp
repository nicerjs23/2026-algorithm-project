// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할2 (김세훈)
//  Trap Read 생성기: 3가지 조건의 read 데이터셋 자동 생성
//
//  [출력 파일 3종]
//   reads_baseline.txt   : 에러 없는 깨끗한 read
//   reads_indel.txt      : InDel 1~2개만 주입 (치환 없음)
//   reads_end_heavy.txt  : 앞 20bp 완벽 + 뒤 10bp에 에러 2개 집중
//
//  [실험 목적]
//   Baseline   → 순수 탐색 속도 비교
//   InDel      → InDel 발생 시 brute-force 정확도 0% 예상, 메인 알고리즘 방어력 검증
//   End-Heavy  → 앞 seed 완벽 보장 → Seed-and-Extend 전략 타당성 증명
//
//  [출력 형식] read_id \t start_pos \t sequence
//   - start_pos : original 기준 0-based 정답 위치
//   - sequence  : InDel 적용 시 길이가 30이 아닐 수 있음
//
//  [팀원 참고]
//   SNP 비율별 reference 4종 (이동건) × read 조건 3종 (김세훈) = 12개 실험
//   SEED=42 고정 → 팀원 누구나 동일한 파일 재현 가능
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <array>

// ──────────────── 파라미터 ────────────────
static const std::string IN_ORIGINAL    = "original_synthetic_1M.txt";

static const std::string OUT_BASELINE   = "reads_baseline.txt";
static const std::string OUT_INDEL      = "reads_indel.txt";
static const std::string OUT_END_HEAVY  = "reads_end_heavy.txt";

static const int      L     = 30;      // read 길이
static const int      M     = 100000;  // read 수
static const int      SEED_BPE = 20;   // Seed 구간 (앞 20bp, 에러 없음)
static const int      TAIL  = L - SEED_BPE;  // Tail 구간 (뒤 10bp)
static const unsigned SEED  = 42;      // 고정 시드
// ──────────────────────────────────────────

static const std::array<char, 4> BASES = {'A', 'C', 'G', 'T'};

static char randomBase(std::mt19937& rng) {
    return BASES[std::uniform_int_distribution<int>(0, 3)(rng)];
}

static char mutate(char c, std::mt19937& rng) {
    int idx = (c=='A')?0:(c=='C')?1:(c=='G')?2:3;
    return BASES[(idx + std::uniform_int_distribution<int>(1,3)(rng)) & 3];
}

// ── Baseline: 에러 없이 그대로 자르기 ─────────────────────────
static void genBaseline(const std::string& orig, std::mt19937& rng) {
    std::ofstream fout(OUT_BASELINE);
    fout << "read_id\tstart_pos\tsequence\n";
    std::uniform_int_distribution<int> pos(0, (int)orig.size() - L);
    for (int i = 0; i < M; i++) {
        int start = pos(rng);
        fout << i << '\t' << start << '\t' << orig.substr(start, L) << '\n';
    }
    std::cout << "[완료] " << OUT_BASELINE << " (에러 없음)\n";
}

// ── InDel Focus: InDel 1~2개 주입, 치환 없음 ──────────────────
// 실험 목적: InDel로 인한 길이 변화가 brute-force를 무력화하는지 검증
static void genInDel(const std::string& orig, std::mt19937& rng) {
    std::ofstream fout(OUT_INDEL);
    fout << "read_id\tstart_pos\tsequence\n";
    std::uniform_int_distribution<int> pos(0, (int)orig.size() - L);
    std::uniform_int_distribution<int> indelCount(1, 2); // InDel 1~2개
    std::uniform_int_distribution<int> indelType(0, 1);  // 0=삽입, 1=삭제

    for (int i = 0; i < M; i++) {
        int start = pos(rng);
        std::string seg = orig.substr(start, L);

        int n = indelCount(rng);
        std::uniform_int_distribution<int> indelPos(0, (int)seg.size() - 1);

        for (int k = 0; k < n; k++) {
            int p = indelPos(rng) % (int)seg.size();
            if (indelType(rng) == 0) {
                // Insertion: 랜덤 염기 삽입
                seg.insert(seg.begin() + p, randomBase(rng));
            } else {
                // Deletion: 해당 위치 염기 삭제
                if (seg.size() > 1) seg.erase(seg.begin() + p);
            }
        }

        fout << i << '\t' << start << '\t' << seg << '\n';
    }
    std::cout << "[완료] " << OUT_INDEL << " (InDel 1~2개, 치환 없음)\n";
}

// ── End-Heavy: 앞 20bp 완벽, 뒤 10bp에 에러 2개 집중 ──────────
// 실험 목적: 앞 seed(20bp) exact match → Seed-and-Extend 전략 타당성 증명
static void genEndHeavy(const std::string& orig, std::mt19937& rng) {
    std::ofstream fout(OUT_END_HEAVY);
    fout << "read_id\tstart_pos\tsequence\n";
    std::uniform_int_distribution<int> pos(0, (int)orig.size() - L);
    std::uniform_int_distribution<int> tailPos(SEED_BPE, L - 1); // 뒤 10bp 내 위치
    std::uniform_int_distribution<int> errType(0, 2); // 0=치환, 1=삽입, 2=삭제

    for (int i = 0; i < M; i++) {
        int start = pos(rng);
        std::string seg = orig.substr(start, L);

        // 뒤 10bp에 에러 정확히 2개 주입 (서로 다른 위치)
        std::vector<int> usedPos;
        int injected = 0;
        int attempts = 0;

        while (injected < 2 && attempts < 20) {
            attempts++;
            int p = tailPos(rng);
            if (p >= (int)seg.size()) continue;
            bool dup = false;
            for (int u : usedPos) if (u == p) { dup = true; break; }
            if (dup) continue;

            int et = errType(rng);
            if (et == 0) {
                // 치환
                seg[p] = mutate(seg[p], rng);
            } else if (et == 1) {
                // 삽입
                seg.insert(seg.begin() + p, randomBase(rng));
            } else {
                // 삭제
                if (seg.size() > SEED_BPE + 1) seg.erase(seg.begin() + p);
                else seg[p] = mutate(seg[p], rng); // 너무 짧아지면 치환으로 대체
            }
            usedPos.push_back(p);
            injected++;
        }

        fout << i << '\t' << start << '\t' << seg << '\n';
    }
    std::cout << "[완료] " << OUT_END_HEAVY << " (앞 20bp 완벽, 뒤 10bp 에러 2개)\n";
}

int main() {
    std::ifstream fin(IN_ORIGINAL);
    if (!fin.is_open()) {
        std::cerr << "[오류] " << IN_ORIGINAL << " 없음.\n";
        return 1;
    }
    std::string orig;
    std::getline(fin, orig);
    fin.close();
    std::cout << "[로드] " << IN_ORIGINAL << " (" << orig.size() << " 염기)\n\n";

    // 각 조건마다 독립적인 시드로 재현 가능하게
    std::mt19937 rng1(SEED);
    std::mt19937 rng2(SEED + 1);
    std::mt19937 rng3(SEED + 2);

    genBaseline(orig,  rng1);
    genInDel(orig,     rng2);
    genEndHeavy(orig,  rng3);

    std::cout << "\n[전체 완료] 3가지 조건 read 생성\n";
    std::cout << "  SNP 비율별 reference 4종 (이동건) × 3종 = 12개 실험 준비 완료\n";
    return 0;
}
