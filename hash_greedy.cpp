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
//  [bruteforce 와의 비교 포인트]
//   - bruteforce: O(N * M * L) = ~3 * 10^12 연산
//   - hash+greedy: O(N) 인덱스 구축 + read 당 평균 O(후보수 * L)
//     -> 후보수가 평균 ~1개 가까이면 사실상 O(M * L)
//
//  [측정 항목]
//   - Index 구축 시간 / Search 시간 / 총 메모리
//   - Mapping rate    : 매핑 성공 read / 전체 read
//   - Position 정확도 : 정답 start_pos 와 일치한 read / 전체 read
//   - 재구축 일치율    : reconstructed vs original 글자 일치율
// =============================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <ctime>

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
// ──────────────────────────────────────────

struct Read {
    int start_pos;   // 정답 위치 (reads.txt 에 기록된 0-based)
    string seq;
};

double check_memory() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
#else
    // POSIX (macOS / Linux): getrusage 로 최대 RSS 측정
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
  #ifdef __APPLE__
    // macOS: ru_maxrss 는 바이트 단위
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
  #else
    // Linux: ru_maxrss 는 킬로바이트 단위
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

// Reference 에 대해 k-mer 해시 인덱스 구축
// key   : 길이 k 의 substring
// value : reference 내에서 해당 substring 이 시작되는 위치들
unordered_map<string, vector<int>> buildHashIndex(const string& ref, int k) {
    unordered_map<string, vector<int>> index;
    index.reserve(ref.size());  // 충돌 줄이기

    const int limit = (int)ref.size() - k;
    for (int i = 0; i <= limit; ++i) {
        index[ref.substr(i, k)].push_back(i);
    }
    return index;
}

// Hash + Greedy 로 read 의 매핑 위치 찾기
// 반환: best 위치 (없으면 -1)
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
    // mismatch <= D 이면 셋 중 적어도 하나는 perfect match -> hash 조회로 잡힘
    const int num_seeds = max_mismatch + 1;
    for (int s = 0; s < num_seeds; ++s) {
        const int seed_offset = s * k;
        if (seed_offset + k > L) break;

        const string seed = read.substr(seed_offset, k);
        auto it = index.find(seed);
        if (it == index.end()) continue;

        // 각 후보 위치에서 read 전체를 비교 (extend 단계)
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

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8
#endif

    cout << "데이터 로딩 중...\n";
    string reference_genome = readReference("reference_synthetic.txt");
    vector<Read> reads      = readReads("reads_synthetic.txt");
    string original_seq     = readReference("original_synthetic_1M.txt");

    //string reference_genome = readReference("reference_yeast.txt");
    //vector<Read> reads      = readReads("reads_yeast.txt");
    //string original_seq     = readReference("original_yeast_1M.txt");

    if (reference_genome.empty() || reads.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    cout << "Hash Table + Greedy mapping...\n";

    // ──────────────── [1/2] 인덱스 구축 ────────────────
    cout << "  [1/2] Hash index 구축 (k=" << SEED_LEN << ") ...\n";
    clock_t t_index_start = clock();
    auto index = buildHashIndex(reference_genome, SEED_LEN);
    clock_t t_index_end = clock();
    double index_time = (double)(t_index_end - t_index_start) / CLOCKS_PER_SEC;
    cout << "        unique " << SEED_LEN << "-mer: " << index.size() << " 개\n";

    // 평균/최대 후보 수 (repeat 영향 측정 지표)
    size_t bucket_max = 0, bucket_sum = 0;
    for (auto& kv : index) {
        bucket_sum += kv.second.size();
        if (kv.second.size() > bucket_max) bucket_max = kv.second.size();
    }
    double bucket_avg = (double)bucket_sum / index.size();

    // ──────────────── [2/2] read 매핑 ────────────────
    cout << "  [2/2] read 매핑 중...\n";
    clock_t t_search_start = clock();

    const int N = (int)reference_genome.size();
    string reconstructed_seq(N, '-');

    int total            = (int)reads.size();
    int mapped           = 0;  // 매핑 성공
    int correct_position = 0;  // 정답 start_pos 와 일치

    const int progress_step = total / 10;
    for (int i = 0; i < total; ++i) {
        int pos = hashGreedySearch(reference_genome, index, reads[i].seq,
                                   SEED_LEN, MAX_MISMATCH);

        if (pos != -1) {
            ++mapped;
            if (pos == reads[i].start_pos) ++correct_position;

            const int L = (int)reads[i].seq.size();
            for (int j = 0; j < L; ++j) {
                if (pos + j < N)
                    reconstructed_seq[pos + j] = reads[i].seq[j];
            }
        }

        if (progress_step > 0 && (i + 1) % progress_step == 0) {
            cout << "  Progress: " << (i + 1) << " / " << total
                 << " (" << (100 * (i + 1) / total) << "%)\n";
        }
    }

    clock_t t_search_end = clock();
    double search_time = (double)(t_search_end - t_search_start) / CLOCKS_PER_SEC;

    // ──────────────── 재구축 일치율 ────────────────
    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original_seq[i])
            ++mismatched;
    }

    double memory       = check_memory();
    double total_time   = index_time + search_time;
    double mapping_rate = 100.0 * mapped           / total;
    double pos_acc      = 100.0 * correct_position / total;
    double correct_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n=========== Hash Table + Greedy 결과 ===========\n";
    cout << "[시간]\n";
    cout << "  Index 구축 시간     : " << index_time   << " 초\n";
    cout << "  Search 시간         : " << search_time  << " 초\n";
    cout << "  총 실행 시간        : " << total_time   << " 초\n";
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
    cout << "================================================\n";

    // 재구축된 염기서열 파일 저장
    ofstream fout("reconstructed_hash_greedy.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_hash_greedy\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed_seq.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_hash_greedy.txt\n";
    } else {
        cerr << "Error: reconstructed_hash_greedy.txt 저장 실패\n";
    }

    return 0;
}
