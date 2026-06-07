// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할3 (이동건)
//  Hash Table + Greedy (Seed-and-Extend, 비둘기집 원리)
//
//  [알고리즘 개요]
//   1) Reference 를 k-mer (k = L/(D+1) = 10) 단위로 해시 인덱싱
//      -> unordered_map<string, vector<int>> : k-mer -> 등장 위치 리스트
//   2) 각 read 를 (D+1)=3 등분 -> 각 조각을 seed 로 hash 조회
//      -> 비둘기집 원리: mismatch <= 2 이면 셋 중 하나는 perfect match
//   3) 후보 위치마다 read 전체를 비교 (extend) 하여 mismatch 카운트
//      -> mismatch 가 가장 작은 위치를 Greedy 하게 선택
//
//  [InDel 에 대한 의도된 약점]
//   - Greedy 는 "글자 단위 1:1 일치 카운트" 만 한다.
//   - read 에 삽입/결실이 생기면 그 뒤 글자가 통째로 밀려 mismatch 가 폭증한다.
//   -> InDel 환경에서 무너지도록 설계 (보고서 비교축: Greedy vs DP).
//
//  [실험 자동화 - 한 번 실행으로 전체 매트릭스]
//   dataset_generator 가 만든 SNP별 reference + 공용 reads 를 사용해
//   SNP 4종(0.1/1/2/5%) × 데이터셋 4종(baseline/indel/end_heavy/yeast)을
//   순서대로 모두 매핑하고 results_hash.csv 에 누적 저장한다.
//   (reference 만 SNP별, reads 는 SNP 무관이라 공용 1벌)
//   같은 reference 를 쓰는 synthetic 3종은 인덱스를 1번만 구축해 재사용.
//
//  [측정 항목 (다른 알고리즘 파일과 동일 포맷)]
//   - 걸린 시간 / 메모리 / N / 불일치 글자 수 / 재구축 일치율(정확도)
// =============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <ctime>
#include <iomanip>

// ---- 플랫폼별 메모리 측정 / 콘솔 인코딩 분기 ----
#ifdef _WIN32
  #include <windows.h>
  #include <psapi.h>
#else
  #include <sys/resource.h>
#endif

using namespace std;

// ──────────────── 파라미터 ────────────────
static const int SEED_LEN     = 10;  // k-mer 길이 (L=30, D=2 -> L/(D+1)=10)
static const int MAX_MISMATCH = 2;   // 허용 mismatch 한계

static const string CSV_FILE   = "results_hash.csv";
static const string ALGO_NAME  = "hash_greedy";
// ──────────────────────────────────────────

struct Read {
    int start_pos;   // 정답 위치 (reads.txt 에 기록된 0-based)
    string seq;
};

// SNP 비율 한 단계 (라벨 + CSV 표기)
struct SnpLevel { string label; string pct; };

// 데이터셋 정의 (reference 는 SNP별, reads/original 은 공용)
struct DatasetDef {
    string label;        // CSV dataset 컬럼
    string ref_prefix;   // reference_<prefix>_snp<label>.txt
    string reads;        // 공용 reads 파일
    string original;     // 정답 원본
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
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);   // macOS: 바이트
  #else
    return (double)ru.ru_maxrss / 1024.0;              // Linux: 킬로바이트
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

// Reference 에 대해 k-mer 해시 인덱스 구축
unordered_map<string, vector<int>> buildHashIndex(const string& ref, int k) {
    unordered_map<string, vector<int>> index;
    index.reserve(ref.size());
    const int limit = (int)ref.size() - k;
    for (int i = 0; i <= limit; ++i) {
        index[ref.substr(i, k)].push_back(i);
    }
    return index;
}

// Hash + Greedy 로 read 의 매핑 위치 찾기 (없으면 -1)
int hashGreedySearch(
    const string& ref,
    const unordered_map<string, vector<int>>& index,
    const string& read,
    int k,
    int max_mismatch)
{
    const int L = (int)read.size();
    const int N = (int)ref.size();

    int best_pos = -1;
    int best_mm  = max_mismatch + 1;

    // 비둘기집 원리: read 를 (D+1) 등분하여 각 조각을 seed 로 사용
    const int num_seeds = max_mismatch + 1;
    for (int s = 0; s < num_seeds; ++s) {
        const int seed_offset = s * k;
        if (seed_offset + k > L) break;

        const string seed = read.substr(seed_offset, k);
        auto it = index.find(seed);
        if (it == index.end()) continue;

        for (int hit : it->second) {
            const int cand = hit - seed_offset;
            if (cand < 0 || cand + L > N) continue;

            int mm = 0;
            for (int j = 0; j < L; ++j) {
                if (ref[cand + j] != read[j]) {
                    if (++mm >= best_mm) break;  // 가지치기
                }
            }

            if (mm < best_mm) {
                best_mm  = mm;
                best_pos = cand;
                if (best_mm == 0) return best_pos;  // perfect -> 조기 종료 (greedy)
            }
        }
    }
    return best_pos;
}

// 결과 1줄을 results_hash.csv 에 누적 (bruteforce.csv 와 동일 스키마)
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
    for (int i = 0; i < total; ++i) {
        int pos = hashGreedySearch(reference, index, reads[i].seq, SEED_LEN, MAX_MISMATCH);
        if (pos != -1) {
            ++mapped;
            if (pos == reads[i].start_pos) ++correct_position;

            // Greedy 재구축: 1:1 선형 복사 (InDel 보정 없음 = 의도된 약점)
            const int L = (int)reads[i].seq.size();
            for (int j = 0; j < L; ++j)
                if (pos + j < N) reconstructed[pos + j] = reads[i].seq[j];
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
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8
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

    cout << "=== Hash Table + Greedy 전체 매트릭스 실험 ===\n";

    for (const auto& snp : snp_levels) {
        // 같은 reference 파일은 인덱스를 1번만 구축해 재사용
        string cur_ref_file;
        string reference;
        unordered_map<string, vector<int>> index;
        double index_time = 0.0;

        for (const auto& ds : datasets) {
            string ref_file = ds.ref_prefix + "_snp" + snp.label + ".txt";

            // reference 캐시 갱신 (파일이 바뀔 때만 로드/인덱싱)
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
                index = buildHashIndex(reference, SEED_LEN);
                index_time = (double)(clock() - t) / CLOCKS_PER_SEC;
                cur_ref_file = ref_file;
            } else if (reference.empty()) {
                continue;  // 직전 reference 로드 실패 상태
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
