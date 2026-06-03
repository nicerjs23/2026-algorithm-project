// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할3 (이동건)
//  Hash Table + DP (Seed-and-Extend, 비둘기집 원리, Needleman-Wunsch 점수)
//
//  [Hash + Greedy 와의 차이]
//   - Greedy 는 "글자 단위 일치 카운트" 만 함 -> InDel 못 잡음 (글자가 밀려서 mismatch 폭증)
//   - DP   는 NW-style 점수로 후보를 평가 -> 삽입/결실 1~2 개까지 흡수 가능
//   - 단, DP 는 후보당 O(L*window) 의 행렬 계산이 들어가서 Greedy 대비 명백히 느림
//   -> 보고서 핵심 비교축: "인델 환경에서 Greedy 무너지고 DP 살아남음"
//
//  [RB-tree + DP (팀원2) 와의 차이]
//   - 인덱스 자료구조만 다름: RB-tree(O(log n) 조회) vs Hash(O(1) 평균 조회)
//   - DP 로직은 동일 -> 정확도는 거의 같고, 인덱스 시간/메모리만 차이
//   -> 보고서 부가 비교축: 자료구조에 따른 인덱스 속도 차이
//
//  [DP 점수 체계]
//   MATCH=+2, MISMATCH=-1, GAP=-2
//   변이 1개당 손실: mismatch -3점 (+2 -> -1), indel -4점 (+2 -> -2)
//   SCORE_THRESHOLD = READ_LEN*MATCH
//                   + MAX_MISMATCH * (MISMATCH - MATCH)
//                   + MAX_INDEL    * (GAP      - MATCH)
//                   = 30*2 + 2*(-3) + 2*(-4) = 60 - 6 - 8 = 46
//   -> "30bp 중 mismatch 2개 + indel 2개 동시 발생까지 통과" 라는 의미
// =============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <ctime>
#include <algorithm>

// ---- 플랫폼별 메모리 측정 / 콘솔 인코딩 분기 ----
#ifdef _WIN32
  #include <windows.h>
  #include <psapi.h>
#else
  #include <sys/resource.h>
#endif

using namespace std;

// ──────────────── 파라미터 ────────────────
static const int K_MER        = 10;  // seed 길이 (L/(D+1) = 30/3)
static const int READ_LEN     = 30;
static const int MAX_MISMATCH = 2;   // 허용 mismatch (reference 측 SNP 대응)
static const int MAX_INDEL    = 2;   // 허용 InDel    (read 측 시퀀서 오류 대응)

// DP 점수 체계
static const int SCORE_MATCH    =  2;
static const int SCORE_MISMATCH = -1;
static const int SCORE_GAP      = -2;

// 통과 임계값: mismatch / indel 둘 다를 비용으로 반영
//   완벽 정렬       = 30*MATCH                                = 60
//   - mismatch 1개  = MATCH 잃고 MISMATCH 로     -> -(MATCH-MISMATCH) = -3
//   - indel 1개     = MATCH 잃고 GAP 로          -> -(MATCH-GAP)      = -4
// SCORE_THRESHOLD = 60 - 2*3 - 2*4 = 60 - 14 = 46
//   "mismatch 2 개 + indel 2 개까지 통과" 라는 의미
static const int SCORE_THRESHOLD = (READ_LEN * SCORE_MATCH)
                                 + MAX_MISMATCH * (SCORE_MISMATCH - SCORE_MATCH)
                                 + MAX_INDEL    * (SCORE_GAP      - SCORE_MATCH);

// DP window 여유 (인델 발생 시 reference 측 길이 보정)
static const int WINDOW_PADDING = 5;
// ──────────────────────────────────────────

struct Read {
    int start_pos;
    string seq;
};

double check_memory() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
#else
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
  #ifdef __APPLE__
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
  #else
    return (double)ru.ru_maxrss / 1024.0;
  #endif
#endif
}

string readReference(const string& filename) {
    ifstream input(filename);
    if (!input.is_open()) {
        cerr << "Error " << filename << " open fail.\n";
        return "";
    }
    string content;
    getline(input, content);
    return content;
}

vector<Read> readReads(const string& filename) {
    ifstream input(filename);
    vector<Read> reads;
    if (!input.is_open()) {
        cerr << "Error " << filename << " open fail.\n";
        return reads;
    }
    reads.reserve(100000);
    string line;
    getline(input, line);  // 헤더 스킵
    while (getline(input, line)) {
        size_t tab1 = line.find('\t');
        size_t tab2 = line.find('\t', tab1 + 1);
        Read r;
        r.start_pos = stoi(line.substr(tab1 + 1, tab2 - tab1 - 1));
        r.seq       = line.substr(tab2 + 1);
        reads.push_back(r);
    }
    return reads;
}

// Reference 의 k-mer 해시 인덱스 (hash_greedy 와 동일 구조)
unordered_map<string, vector<int>> buildHashIndex(const string& ref, int k) {
    unordered_map<string, vector<int>> index;
    index.reserve(ref.size());
    const int limit = (int)ref.size() - k;
    for (int i = 0; i <= limit; ++i) {
        index[ref.substr(i, k)].push_back(i);
    }
    return index;
}

// Needleman-Wunsch 스타일 DP (1D 배열 메모리 최적화).
// read 와 ref_window 정렬 점수의 최댓값을 반환한다.
int calculateDPScore(const string& read_seq, const string& ref_window) {
    const int n = (int)read_seq.size();
    const int m = (int)ref_window.size();

    vector<int> prev_row(m + 1, 0);
    vector<int> curr_row(m + 1, 0);

    for (int j = 1; j <= m; ++j) prev_row[j] = j * SCORE_GAP;

    for (int i = 1; i <= n; ++i) {
        curr_row[0] = i * SCORE_GAP;
        for (int j = 1; j <= m; ++j) {
            int diag = prev_row[j - 1]
                     + (read_seq[i - 1] == ref_window[j - 1] ? SCORE_MATCH : SCORE_MISMATCH);
            int up   = prev_row[j]     + SCORE_GAP;
            int left = curr_row[j - 1] + SCORE_GAP;
            curr_row[j] = max({diag, up, left});
        }
        prev_row = curr_row;
    }

    int best = prev_row[0];
    for (int j = 1; j <= m; ++j) best = max(best, prev_row[j]);
    return best;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif

    // ────── 데이터셋 선택 (사용할 한 묶음만 주석 풀고 나머지는 닫기) ──────
    // [기본] 인공 서열
    string reference_file = "reference_synthetic.txt";
    string reads_file     = "reads_synthetic.txt";
    string original_file  = "original_synthetic_1M.txt";

    // [비교] 빵효모 (반복서열 영향 확인)
    //string reference_file = "reference_yeast.txt";
    //string reads_file     = "reads_yeast.txt";
    //string original_file  = "original_yeast_1M.txt";

    // [Trap] Baseline (에러 없는 깨끗한 read)
    //string reference_file = "reference_synthetic.txt";
    //string reads_file     = "reads_baseline.txt";
    //string original_file  = "original_synthetic_1M.txt";

    // [Trap] InDel (삽입/결실 폭격 - DP 가 빛나는 환경)
    //string reference_file = "reference_synthetic.txt";
    //string reads_file     = "reads_indel.txt";
    //string original_file  = "original_synthetic_1M.txt";

    // [Trap] End-Heavy (read 후반부에 에러 집중)
    //string reference_file = "reference_synthetic.txt";
    //string reads_file     = "reads_end_heavy.txt";
    //string original_file  = "original_synthetic_1M.txt";
    // ────────────────────────────────────────────────────────────

    cout << "데이터 로딩 중...\n";
    cout << "  Reference: " << reference_file << "\n";
    cout << "  Reads    : " << reads_file     << "\n";
    string reference_genome = readReference(reference_file);
    vector<Read> reads      = readReads(reads_file);
    string original_seq     = readReference(original_file);

    if (reference_genome.empty() || reads.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    cout << "Hash Table + DP mapping...\n";

    // ──────────────── [1/2] 인덱스 구축 ────────────────
    cout << "  [1/2] Hash index 구축 (k=" << K_MER << ") ...\n";
    clock_t t_idx_s = clock();
    auto index = buildHashIndex(reference_genome, K_MER);
    clock_t t_idx_e = clock();
    double idx_time = (double)(t_idx_e - t_idx_s) / CLOCKS_PER_SEC;
    cout << "        unique " << K_MER << "-mer: " << index.size() << " 개\n";

    size_t bucket_sum = 0, bucket_max = 0;
    for (auto& kv : index) {
        bucket_sum += kv.second.size();
        if (kv.second.size() > bucket_max) bucket_max = kv.second.size();
    }
    double bucket_avg = (double)bucket_sum / index.size();

    // ──────────────── [2/2] read 매핑 (DP 점수 평가) ────────────────
    cout << "  [2/2] read 매핑 중 (DP 점수 평가) ...\n";
    clock_t t_srch_s = clock();

    const int N = (int)reference_genome.size();
    string reconstructed_seq(N, '-');

    int total            = (int)reads.size();
    int mapped           = 0;
    int correct_position = 0;

    const int progress_step = total / 10;
    for (int i = 0; i < total; ++i) {
        int best_pos   = -1;
        int best_score = SCORE_THRESHOLD - 1;  // 임계값 미만이면 매핑 실패

        const int L = (int)reads[i].seq.size();

        // 비둘기집 원리: read 를 (D+1)=3 등분, 각 조각을 seed 로 사용
        for (int offset = 0; offset + K_MER <= L && offset <= MAX_MISMATCH * K_MER; offset += K_MER) {
            const string seed = reads[i].seq.substr(offset, K_MER);
            auto it = index.find(seed);
            if (it == index.end()) continue;

            for (int pos : it->second) {
                const int ref_start = pos - offset;
                if (ref_start < 0 || ref_start >= N) continue;

                // 인델 흡수용 window: read 길이보다 약간 길게
                const int window_len = min(READ_LEN + WINDOW_PADDING, N - ref_start);
                const string ref_window = reference_genome.substr(ref_start, window_len);

                const int score = calculateDPScore(reads[i].seq, ref_window);

                if (score > best_score) {
                    best_score = score;
                    best_pos   = ref_start;
                }
            }
        }

        if (best_pos != -1) {
            ++mapped;
            if (best_pos == reads[i].start_pos) ++correct_position;

            for (int j = 0; j < L; ++j) {
                if (best_pos + j < N)
                    reconstructed_seq[best_pos + j] = reads[i].seq[j];
            }
        }

        if (progress_step > 0 && (i + 1) % progress_step == 0) {
            cout << "  Progress: " << (i + 1) << " / " << total
                 << " (" << (100 * (i + 1) / total) << "%)\n";
        }
    }

    clock_t t_srch_e = clock();
    double srch_time = (double)(t_srch_e - t_srch_s) / CLOCKS_PER_SEC;

    // ──────────────── 재구축 일치율 ────────────────
    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original_seq[i])
            ++mismatched;
    }

    double memory       = check_memory();
    double total_time   = idx_time + srch_time;
    double mapping_rate = 100.0 * mapped           / total;
    double pos_acc      = 100.0 * correct_position / total;
    double correct_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n=========== Hash Table + DP 결과 ===========\n";
    cout << "[데이터셋]\n";
    cout << "  Reference : " << reference_file << "\n";
    cout << "  Reads     : " << reads_file     << "\n";
    cout << "[시간]\n";
    cout << "  Index 구축 시간     : " << idx_time   << " 초\n";
    cout << "  Search (DP) 시간    : " << srch_time  << " 초\n";
    cout << "  총 실행 시간        : " << total_time << " 초\n";
    cout << "[메모리]\n";
    cout << "  사용 중 메모리      : " << memory       << " MB\n";
    cout << "[인덱스 통계]\n";
    cout << "  unique k-mer 수    : " << index.size() << "\n";
    cout << "  k-mer 당 평균 위치  : " << bucket_avg   << " 개\n";
    cout << "  k-mer 당 최대 위치  : " << bucket_max   << " 개 (repeat 지표)\n";
    cout << "[정확도]\n";
    cout << "  총 read 수         : " << total            << "\n";
    cout << "  매핑 성공 read     : " << mapped           << " (" << mapping_rate << "%)\n";
    cout << "  정답 위치 일치     : " << correct_position << " (" << pos_acc      << "%)\n";
    cout << "  재구축 일치율      : " << correct_rate     << "%\n";
    cout << "  불일치 글자 수     : " << mismatched       << " / " << N << "\n";
    cout << "[DP 설정]\n";
    cout << "  점수: MATCH=+" << SCORE_MATCH << " MISMATCH=" << SCORE_MISMATCH
         << " GAP=" << SCORE_GAP << "\n";
    cout << "  SCORE_THRESHOLD    : " << SCORE_THRESHOLD
         << " (mismatch " << MAX_MISMATCH << "개 + indel " << MAX_INDEL << "개 허용 기준)\n";
    cout << "============================================\n";

    // 재구축된 염기서열 파일 저장
    ofstream fout("reconstructed_hash_dp.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_hash_dp\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed_seq.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_hash_dp.txt\n";
    } else {
        cerr << "Error: reconstructed_hash_dp.txt 저장 실패\n";
    }

    return 0;
}
