// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할2 (김세훈)
//  BWT + FM-index 기반 Read 매핑
//
//  [알고리즘 흐름]
//   1단계: Suffix Array 구축 (prefix doubling, O(N log²N))
//   2단계: BWT 생성 (SA로 각 접미사 앞 글자 추출)
//   3단계: FM-index 구축 (C 테이블 + Occ 테이블)
//   4단계: Read 매핑 (backward search)
//   5단계: Mismatch 처리 (비둘기집 원리, D=2 → 3구간 각 10bp)
//   6단계: 정확도 측정 (found_pos vs true_pos)
//
//  [비둘기집 원리]
//   L=30, D=2 → 3구간으로 분할 (각 10bp)
//   최소 1구간은 반드시 exact match → 그 구간으로 후보 위치 탐색
//   → 후보 위치에서 전체 read mismatch ≤ D 확인
//
//  [입력]
//   reference_genome.txt : SNP 심긴 레퍼런스 (1번 이동건)
//   reads.txt            : read_id / start_pos / sequence (2번 김세훈)
//
//  [출력]
//   result_bwt.txt : read_id / true_pos / found_pos / mismatch / status
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <chrono>

// ──────────────── 파라미터 ────────────────
static const std::string IN_REFERENCE = "reference_genome.txt";
static const std::string IN_READS     = "reads.txt";
static const std::string OUT_RESULT   = "result_bwt.txt";

static const int D        = 2;    // 최대 허용 mismatch
static const int L        = 30;   // read 길이
static const int SEG_LEN  = L / (D + 1);   // 비둘기집 구간 길이 = 10
static const int MAX_CAND = 200;  // 구간당 최대 후보 수 (반복 서열 폭주 방지)
static const int ALPHA    = 5;    // 알파벳 크기 ($=0, A=1, C=2, G=3, T=4)
// ─────────────────────────────────────────

// 문자 → 인덱스
static int cidx(char c) {
    switch (c) {
        case '$': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        default:  return -1;
    }
}

// ── 1단계: Suffix Array (prefix doubling) ──────────────────────
// 쌍(rank[i], rank[i+gap])으로 정렬을 반복해 O(N log²N) 구축
std::vector<int> buildSA(const std::string& s) {
    int n = (int)s.size();
    std::vector<int> sa(n), rank_(n), tmp(n);
    std::iota(sa.begin(), sa.end(), 0);
    for (int i = 0; i < n; i++) rank_[i] = (unsigned char)s[i];

    for (int gap = 1; gap < n; gap <<= 1) {
        // (rank[i], rank[i+gap]) 쌍으로 정렬
        std::sort(sa.begin(), sa.end(), [&](int a, int b) {
            if (rank_[a] != rank_[b]) return rank_[a] < rank_[b];
            int ra = (a + gap < n) ? rank_[a + gap] : -1;
            int rb = (b + gap < n) ? rank_[b + gap] : -1;
            return ra < rb;
        });
        // 새 rank 계산
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            tmp[sa[i]] = tmp[sa[i - 1]];
            int ra1 = rank_[sa[i]],     rb1 = rank_[sa[i - 1]];
            int ra2 = (sa[i]     + gap < n) ? rank_[sa[i]     + gap] : -1;
            int rb2 = (sa[i - 1] + gap < n) ? rank_[sa[i - 1] + gap] : -1;
            if (ra1 != rb1 || ra2 != rb2) tmp[sa[i]]++;
        }
        rank_ = tmp;
        if (rank_[sa[n - 1]] == n - 1) break;
    }
    return sa;
}

// ── 2단계: BWT 생성 ────────────────────────────────────────────
// SA[i]가 가리키는 접미사의 바로 앞 글자 = BWT[i]
std::string buildBWT(const std::string& s, const std::vector<int>& sa) {
    int n = (int)s.size();
    std::string bwt(n, ' ');
    for (int i = 0; i < n; i++)
        bwt[i] = (sa[i] == 0) ? '$' : s[sa[i] - 1];
    return bwt;
}

// ── 3단계: FM-index ────────────────────────────────────────────
struct FMIndex {
    std::vector<int>                  C;    // C[ALPHA]: 각 문자보다 작은 문자 수
    std::vector<std::array<int, ALPHA>> Occ; // Occ[i][c]: BWT[0..i-1]에서 c 등장 횟수

    explicit FMIndex(const std::string& bwt) {
        int n = (int)bwt.size();
        // Occ 구축 (prefix count)
        Occ.assign(n + 1, {});
        for (int i = 0; i < n; i++) {
            Occ[i + 1] = Occ[i];
            int c = cidx(bwt[i]);
            if (c >= 0) Occ[i + 1][c]++;
        }
        // C 구축: C[c] = sum of Occ[n][0..c-1]
        C.resize(ALPHA, 0);
        for (int c = 1; c < ALPHA; c++)
            C[c] = C[c - 1] + Occ[n][c - 1];
    }
};

// ── 4단계: Backward Search ─────────────────────────────────────
// 패턴의 정확 매칭 SA 범위 [top, bottom] 반환. top > bottom이면 미매칭.
std::pair<int, int> backwardSearch(const std::string& pattern,
                                   const FMIndex& fm) {
    int n   = (int)fm.Occ.size() - 1;
    int top = 0, bottom = n - 1;

    for (int i = (int)pattern.size() - 1; i >= 0; i--) {
        int c = cidx(pattern[i]);
        if (c < 0 || c == 0) return {1, 0}; // '$' 또는 알 수 없는 문자
        top    = fm.C[c] + fm.Occ[top][c];
        bottom = fm.C[c] + fm.Occ[bottom + 1][c] - 1;
        if (top > bottom) return {1, 0};
    }
    return {top, bottom};
}

// ── 5단계 보조: mismatch 수 계산 ───────────────────────────────
int countMM(const std::string& ref, int pos, const std::string& read) {
    int mm = 0;
    for (int i = 0; i < L; i++) {
        if (ref[pos + i] != read[i]) mm++;
        if (mm > D) return mm; // 조기 종료
    }
    return mm;
}

// ── main ───────────────────────────────────────────────────────
int main() {
    auto t0 = std::chrono::steady_clock::now();

    // 1. reference 로드
    std::ifstream fref(IN_REFERENCE);
    if (!fref.is_open()) {
        std::cerr << "[오류] " << IN_REFERENCE << " 없음.\n"; return 1;
    }
    std::string ref;
    std::getline(fref, ref);
    fref.close();
    std::cout << "[로드] reference 길이: " << ref.size() << " 염기\n";

    // 2. SA + BWT (sentinel '$' 추가)
    std::cout << "[1단계] Suffix Array 구축 중... (수십 초 소요)\n";
    std::string refS = ref + '$';
    std::vector<int> sa = buildSA(refS);

    std::cout << "[2단계] BWT 생성 중...\n";
    std::string bwt = buildBWT(refS, sa);

    // 3. FM-index
    std::cout << "[3단계] FM-index 구축 중...\n";
    FMIndex fm(bwt);

    auto t1 = std::chrono::steady_clock::now();
    std::cout << "        인덱스 완료 ("
              << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
              << "초)\n\n";

    // 4~5. Read 매핑
    std::cout << "[4~5단계] Read 매핑 중...\n";

    std::ifstream freads(IN_READS);
    if (!freads.is_open()) {
        std::cerr << "[오류] " << IN_READS << " 없음.\n"; return 1;
    }
    std::ofstream fout(OUT_RESULT);
    fout << "read_id\ttrue_pos\tfound_pos\tmismatch\tstatus\n";

    std::string line;
    std::getline(freads, line); // 헤더 스킵

    int total = 0, mapped = 0, correct = 0;
    int refLen = (int)ref.size();

    while (std::getline(freads, line)) {
        std::istringstream iss(line);
        int readId, truePos;
        std::string seq;
        iss >> readId >> truePos >> seq;
        if ((int)seq.size() != L) continue;

        total++;
        int foundPos = -1;
        int foundMM  = D + 1;

        // 비둘기집 원리: D+1=3 구간, 각 SEG_LEN=10 bp
        for (int seg = 0; seg <= D && foundPos == -1; seg++) {
            int segStart = seg * SEG_LEN;
            std::string seed = seq.substr(segStart, SEG_LEN);

            auto [top, bottom] = backwardSearch(seed, fm);
            if (top > bottom) continue;

            // 후보 수 제한 (반복 서열 과다 방지)
            int limit = std::min(bottom, top + MAX_CAND - 1);
            for (int k = top; k <= limit; k++) {
                int pos = sa[k] - segStart; // 전체 read 시작 위치 후보
                if (pos < 0 || pos + L > refLen) continue;

                int mm = countMM(ref, pos, seq);
                if (mm <= D && mm < foundMM) {
                    foundMM  = mm;
                    foundPos = pos;
                }
            }
        }

        // 6. 결과 기록
        std::string status;
        if (foundPos < 0) {
            status  = "unmapped";
            foundMM = -1;
        } else {
            mapped++;
            status = (foundPos == truePos) ? "correct" : "wrong";
            if (foundPos == truePos) correct++;
        }

        fout << readId << '\t' << truePos << '\t' << foundPos
             << '\t' << foundMM << '\t' << status << '\n';
    }

    auto t2 = std::chrono::steady_clock::now();

    // 정확도 출력
    std::cout << "\n[결과]\n";
    std::cout << "  전체 read      : " << total   << "\n";
    std::cout << "  매핑 성공      : " << mapped  << " (" << 100.0 * mapped  / total << "%)\n";
    std::cout << "  정확 위치 매핑 : " << correct << " (" << 100.0 * correct / total << "%)\n";
    std::cout << "  총 소요 시간   : "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t0).count()
              << " ms\n";
    std::cout << "  결과 파일      : " << OUT_RESULT << "\n";

    return 0;
}
