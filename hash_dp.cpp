// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할3 (이동건)
//  Hash Table + DP (Seed-and-Extend, 비둘기집 원리, Needleman-Wunsch 점수)
//
//  [Hash + Greedy 와의 차이]
//   - Greedy 는 "글자 단위 일치 카운트" 만 함 -> InDel 못 잡음 (글자가 밀려서 mismatch 폭증)
//   - DP   는 NW-style 점수로 후보를 평가 -> 삽입/결실 1~2 개까지 흡수 가능
//   -> 보고서 핵심 비교축: "인델 환경에서 Greedy 무너지고 DP 살아남음"
//
//  [InDel 정확도 개선 - traceback 기반 재구축]
//   기존 DP 는 점수만 보고 위치를 정한 뒤 reconstructed[pos+j]=read[j] 로 "1:1 선형 복사"를 했다.
//   -> InDel 이 있으면 삽입/결실 이후 글자가 한 칸씩 밀려, 매핑은 성공해도 재구축이 어긋났다.
//   개선: 최적 후보에 대해 DP 행렬을 traceback 해 정렬을 복원하고, 각 read 글자를
//        "정렬된 reference 좌표"에 정확히 배치한다 (삽입 글자는 버리고 결실 구간은 건너뜀).
//        baseline/end_heavy(치환)에서는 gap 이 없어 기존과 동일하게 동작한다.
//
//  [정렬 방식] Fitting alignment (read 를 reference window 의 최적 구간에 끼워맞춤)
//   - 양끝 reference gap 무비용 (free flank) -> seed 주변 어디서 시작/끝나도 OK
//   - read 글자는 모두 소비 (남는 read 글자는 mismatch/insertion 비용으로 반영)
//
//  [DP 점수 체계] MATCH=+2, MISMATCH=-1, GAP=-2
//   SCORE_THRESHOLD = 30*2 + 2*(-3) + 2*(-4) = 46
//   -> "30bp 중 mismatch 2개 + indel 2개까지 통과" 라는 의미
//
//  [실험 자동화 - 한 번 실행으로 전체 매트릭스]
//   dataset_generator 가 만든 SNP별 reference + 공용 reads 로
//   SNP 4종 × 데이터셋 4종(baseline/indel/end_heavy/yeast)을 모두 매핑하고
//   results_hash.csv 에 누적 (hash_greedy 와 동일 스키마, 같은 파일에 합쳐짐).
// =============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <ctime>
#include <iomanip>
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

// 통과 임계값: mismatch 2개 + indel 2개까지 통과 (= 60 - 2*3 - 2*4 = 46)
static const int SCORE_THRESHOLD = (READ_LEN * SCORE_MATCH)
                                 + MAX_MISMATCH * (SCORE_MISMATCH - SCORE_MATCH)
                                 + MAX_INDEL    * (SCORE_GAP      - SCORE_MATCH);

// DP window 여유 (인델로 인한 좌우 길이 변화 흡수: read 양쪽으로 PADDING 만큼)
static const int WINDOW_PADDING = 5;

static const string CSV_FILE   = "results_hash.csv";
static const string ALGO_NAME  = "hash_dp";
// ──────────────────────────────────────────

struct Read {
    int start_pos;
    string seq;
};

struct SnpLevel { string label; string pct; };

struct DatasetDef {
    string label;
    string ref_prefix;
    string reads;
    string original;
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
    if (!input.is_open()) return "";
    string content;
    getline(input, content);
    return content;
}

vector<Read> readReads(const string& filename) {
    ifstream input(filename);
    vector<Read> reads;
    if (!input.is_open()) return reads;
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

unordered_map<string, vector<int>> buildHashIndex(const string& ref, int k) {
    unordered_map<string, vector<int>> index;
    index.reserve(ref.size());
    const int limit = (int)ref.size() - k;
    for (int i = 0; i <= limit; ++i) index[ref.substr(i, k)].push_back(i);
    return index;
}

// ── Fitting alignment 점수만 계산 (1D 메모리 최적화, 후보 선별용) ──
//   윗줄 초기값 0 : 앞쪽 reference gap 무비용 / curr_row[0]=i*GAP : 앞선 read 는 insertion
//   반환 = 마지막 행 최댓값 : 뒤쪽 reference gap 무비용
int calculateDPScore(const string& read_seq, const string& ref_window) {
    const int n = (int)read_seq.size();
    const int m = (int)ref_window.size();

    vector<int> prev_row(m + 1, 0);
    vector<int> curr_row(m + 1, 0);

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

// ── 최적 후보 1개에 대해 full DP + traceback 으로 정렬 복원 ──
//   reconstructed 의 절대 좌표(window_start 기준)에 read 글자를 배치.
//   반환값: read[0] 가 정렬된 reference 절대 좌표 (정답 위치 비교용)
int reconstructByTraceback(const string& read_seq, const string& reference,
                           int window_start, int window_len, string& reconstructed) {
    const int n = (int)read_seq.size();
    const int m = window_len;
    const int N = (int)reference.size();

    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    for (int i = 1; i <= n; ++i) dp[i][0] = i * SCORE_GAP;   // 앞 read = insertion

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            char rc = reference[window_start + j - 1];
            int diag = dp[i - 1][j - 1] + (read_seq[i - 1] == rc ? SCORE_MATCH : SCORE_MISMATCH);
            int up   = dp[i - 1][j]     + SCORE_GAP;
            int left = dp[i][j - 1]     + SCORE_GAP;
            dp[i][j] = max({diag, up, left});
        }
    }

    // 뒤쪽 reference gap 무비용 -> 마지막 행 최댓값 위치를 끝점으로
    int end_j = 0, best = dp[n][0];
    for (int j = 1; j <= m; ++j)
        if (dp[n][j] >= best) { best = dp[n][j]; end_j = j; }

    int i = n, j = end_j;
    int read0_ref = -1;
    while (i > 0) {
        if (j > 0) {
            char rc = reference[window_start + j - 1];
            int diag = dp[i - 1][j - 1] + (read_seq[i - 1] == rc ? SCORE_MATCH : SCORE_MISMATCH);
            if (dp[i][j] == diag) {
                int abs_pos = window_start + j - 1;
                if (abs_pos >= 0 && abs_pos < N) reconstructed[abs_pos] = read_seq[i - 1];
                if (i == 1) read0_ref = abs_pos;
                --i; --j;
                continue;
            }
            if (dp[i][j] == dp[i - 1][j] + SCORE_GAP) {  // insertion -> 버림
                if (i == 1) read0_ref = window_start + j;
                --i;
                continue;
            }
            --j;  // deletion -> 건너뜀
        } else {
            if (i == 1) read0_ref = window_start;
            --i;
        }
    }
    return read0_ref;
}

void appendResultCsv(const string& algorithm, const string& dataset, const string& snp,
                     double total_sec, double memory, double mapped_pct, double recon_pct) {
    bool is_new = false;
    {
        ifstream check(CSV_FILE);
        if (!check.is_open()) is_new = true;
    }
    ofstream csv(CSV_FILE, ios::app);
    if (!csv.is_open()) { cerr << "[경고] " << CSV_FILE << " 저장 실패\n"; return; }
    if (is_new)
        csv << "algorithm,dataset,snp_rate,total_sec,memory_mb,mapped_pct,reconstruct_pct\n";
    csv << fixed << setprecision(2)
        << algorithm << "," << dataset << "," << snp << ","
        << total_sec << "," << memory << "," << mapped_pct << "," << recon_pct << "\n";
}

// 이미 로드된 reference/index + reads/original 로 한 셀(SNP×데이터셋)을 매핑·기록
void mapCell(const string& reference,
             const unordered_map<string, vector<int>>& index,
             double index_time,
             const vector<Read>& reads, const string& original,
             const string& dataset_label, const string& snp_pct) {

    const int N = (int)reference.size();
    string reconstructed(N, '-');

    const int total = (int)reads.size();
    int mapped = 0, correct_position = 0;

    clock_t t_s = clock();
    for (int r = 0; r < total; ++r) {
        const string& read = reads[r].seq;
        const int L = (int)read.size();

        int best_score  = SCORE_THRESHOLD - 1;  // 임계값 미만이면 매핑 실패
        int best_wstart = -1, best_wlen = 0;

        // 비둘기집 원리: read 를 (D+1)=3 등분, 각 조각을 seed 로 사용
        for (int offset = 0; offset + K_MER <= L && offset <= MAX_MISMATCH * K_MER; offset += K_MER) {
            const string seed = read.substr(offset, K_MER);
            auto it = index.find(seed);
            if (it == index.end()) continue;

            for (int pos : it->second) {
                int anchor = pos - offset;
                int wstart = max(0, anchor - WINDOW_PADDING);
                int wlen   = min(READ_LEN + 2 * WINDOW_PADDING, N - wstart);
                if (wlen <= 0) continue;

                const string ref_window = reference.substr(wstart, wlen);
                int score = calculateDPScore(read, ref_window);
                if (score > best_score) {
                    best_score  = score;
                    best_wstart = wstart;
                    best_wlen   = wlen;
                }
            }
        }

        if (best_wstart != -1) {
            ++mapped;
            int read0_ref = reconstructByTraceback(read, reference, best_wstart, best_wlen,
                                                   reconstructed);
            if (read0_ref == reads[r].start_pos) ++correct_position;
        }
    }
    double search_time = (double)(clock() - t_s) / CLOCKS_PER_SEC;

    int mismatched = 0;
    for (int i = 0; i < N; ++i)
        if (reconstructed[i] != original[i]) ++mismatched;

    double memory              = check_memory();
    double elapsed_sec         = index_time + search_time;
    double mapping_rate        = 100.0 * mapped           / total;
    double pos_acc             = 100.0 * correct_position / total;
    double reconstruction_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);
    cout << "\n----- [" << dataset_label << "] SNP " << snp_pct << " -----\n";
    cout << "걸린 시간              : " << elapsed_sec         << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory              << " MB\n";
    cout << "총 원본 염기서열 길이(N): " << N                   << "\n";
    cout << "일치하지 않는 글자 수   : " << mismatched          << " 개 (미복구 빈칸 포함)\n";
    cout << "재구축 일치율(정확도)   : " << reconstruction_rate << "%\n";
    cout << "  (참고) 매핑 성공률    : " << mapping_rate << "%  /  정답 위치 일치율 : " << pos_acc << "%\n";

    appendResultCsv(ALGO_NAME, dataset_label, snp_pct,
                    elapsed_sec, memory, mapping_rate, reconstruction_rate);
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif

    const vector<SnpLevel> snp_levels = {
        {"1", "0.1%"}, {"10", "1.0%"}, {"20", "2.0%"}, {"50", "5.0%"},
    };
    const vector<DatasetDef> datasets = {
        {"baseline",  "reference_synthetic", "reads_baseline.txt",  "original_synthetic_1M.txt"},
        {"indel",     "reference_synthetic", "reads_indel.txt",     "original_synthetic_1M.txt"},
        {"end_heavy", "reference_synthetic", "reads_end_heavy.txt", "original_synthetic_1M.txt"},
        {"yeast",     "reference_yeast",     "reads_yeast.txt",     "original_yeast_1M.txt"},
    };

    cout << "=== Hash Table + DP 전체 매트릭스 실험 ===\n";
    cout << "  SCORE_THRESHOLD = " << SCORE_THRESHOLD
         << " (mismatch " << MAX_MISMATCH << " + indel " << MAX_INDEL << " 허용)\n";

    for (const auto& snp : snp_levels) {
        string cur_ref_file;
        string reference;
        unordered_map<string, vector<int>> index;
        double index_time = 0.0;

        for (const auto& ds : datasets) {
            string ref_file = ds.ref_prefix + "_snp" + snp.label + ".txt";

            if (ref_file != cur_ref_file) {
                string r = readReference(ref_file);
                if (r.empty()) {
                    cout << "\n[건너뜀] " << ref_file << " 없음 -> ["
                         << ds.label << "] SNP " << snp.pct << "\n";
                    cur_ref_file.clear();
                    continue;
                }
                reference = std::move(r);
                clock_t t = clock();
                index = buildHashIndex(reference, K_MER);
                index_time = (double)(clock() - t) / CLOCKS_PER_SEC;
                cur_ref_file = ref_file;
            } else if (reference.empty()) {
                continue;
            }

            vector<Read> reads = readReads(ds.reads);
            string original    = readReference(ds.original);
            if (reads.empty() || original.empty()) {
                cout << "\n[건너뜀] reads/original 없음 -> [" << ds.label
                     << "] SNP " << snp.pct << " (" << ds.reads << " / " << ds.original << ")\n";
                continue;
            }

            mapCell(reference, index, index_time, reads, original, ds.label, snp.pct);
        }
    }

    cout << "\n[완료] 결과가 " << CSV_FILE << " 에 누적 저장되었습니다.\n";
    return 0;
}
