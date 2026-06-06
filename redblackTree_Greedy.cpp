#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <map>
#include <algorithm>
#include <iomanip>

// 메모리 측정
#ifdef _WIN32
  #include <windows.h>
  #include <psapi.h>
#else
  #include <sys/resource.h>
#endif

using namespace std;

const int k = 10;
const int readLen = 30;
const int maxMismatch = 2;  // 허용 mismatch 개수

struct Read {
    int pos;
    string seq;
};

// 메모리 측정 (실제 프로세스 메모리)
double checkMemory() {
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

// 레퍼런스 및 원본 파일 읽기
string readReference(const string& filename) {
    ifstream input(filename);
    string content;
    if (input.is_open()) {
        getline(input, content);
    }
    return content;
}

// Reads 파일 읽기
vector<Read> readReads(const string& filename) {
    ifstream input(filename);
    vector<Read> reads;
    if (!input.is_open()) {
        return reads;
    }
    reads.reserve(100000);
    string line;
    getline(input, line);
    while (getline(input, line)) {
        size_t t1 = line.find('\t');
        size_t t2 = line.find('\t', t1 + 1);
        Read r;
        r.pos = stoi(line.substr(t1 + 1, t2 - t1 - 1));
        r.seq = line.substr(t2 + 1);
        reads.push_back(r);
    }
    return reads;
}

// Red-Black Tree 기반 Seed-and-Extend 구성
map<string, vector<int>> consSeedIdx(const string& reference, int k) {
    map<string, vector<int>> idx;
    int N = (int)reference.length();
    for (int i = 0; i <= N - k; ++i) {
        idx[reference.substr(i, k)].push_back(i);
    }
    return idx;
}

// Greedy 계산 (mismatch 개수)
int countMismatch(const string& ref, int refStart, const string& readSeq) {
    int mm = 0;
    for (int j = 0; j < readLen; ++j) {
        if (ref[refStart + j] != readSeq[j]) {
            ++mm;
            if (mm > maxMismatch) {
                return mm;
            }
        }
    }
    return mm;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8 (한글 깨짐 방지)
#endif

    // =========================================================================
    // [사용자 설정 파트] 실행 전에 현재 돌리는 데이터셋 정보 입력
    // =========================================================================
    string current_algorithm = "redblackTree_Greedy";
    string current_snp       = "1.0%";       // 예: 0.1%, 1.0%, 2.0%, 5.0%
    string current_dataset   = "yeast";      // 예: yeast, baseline, indel, end_heavy
    // =========================================================================

    cout << "데이터 로딩 중...\n";

    // 인공서열
    //string reference_genome = readReference("reference_synthetic.txt");
    //vector<Read> reads      = readReads("reads_baseline.txt");
    //string original_seq     = readReference("original_synthetic_1M.txt");

    //string reference_genome = readReference("reference_synthetic.txt");
    //vector<Read> reads = readReads("reads_indel.txt");
    //string original_seq     = readReference("original_synthetic_1M.txt");

    //string reference_genome = readReference("reference_synthetic.txt");
    //vector<Read> reads = readReads("reads_end_heavy.txt");
    //string original_seq     = readReference("original_synthetic_1M.txt");

    // 효모
    string reference_genome = readReference("reference_yeast.txt");
    vector<Read> reads      = readReads("reads_yeast.txt");
    string original_seq     = readReference("original_yeast_1M.txt");

    if (reference_genome.empty() || reads.empty() || original_seq.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    int N = reference_genome.length();
    string reconstructed(N, '-');

    cout << "Red-Black Tree Seed-and-Extend...\n";
    map<string, vector<int>> seedIdx = consSeedIdx(reference_genome, k);

    cout << "Seed-and-Extend 매핑 (Greedy)...\n";
    clock_t sTime = clock();

    int mappingCount = 0;
    int exactCount = 0;

    for (size_t i = 0; i < reads.size(); ++i) {
        int bestPos = -1;

        // 비둘기집 원리: 0, 10, 20 세 구간에서 seed 추출
        for (int offset = 0; offset <= 20 && bestPos == -1; offset += 10) {
            string seed = reads[i].seq.substr(offset, k);
            auto it = seedIdx.find(seed); // RB-tree에서 mismatch가 없는 구간 찾기

            if (it != seedIdx.end()) {
                for (int pos : it->second) {
                    int refStart = pos - offset;  // seed 위치 역산
                    if (refStart < 0 || refStart + readLen > N) continue;

                    // Greedy: mismatch <= 2이면 첫 번째 성공 위치 바로 선택
                    int mm = countMismatch(reference_genome, refStart, reads[i].seq);
                    if (mm <= maxMismatch) {
                        bestPos = refStart;
                        break;
                    }
                }
            }
        }

        if (bestPos != -1) {
            mappingCount++;
            if (bestPos == reads[i].pos) {
                exactCount++;
            }
            for (int j = 0; j < readLen; ++j) {
                if (bestPos + j < N) {
                    reconstructed[bestPos + j] = reads[i].seq[j];
                }
            }
        }
    }

    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed[i] != original_seq[i]) {
            ++mismatched;
        }
    }

    clock_t end_time = clock();
    double elapsed_sec = (double)(end_time - sTime) / CLOCKS_PER_SEC;
    double memory = checkMemory();
    double reconstruction_rate = 100.0 * (N - mismatched) / N;
    double mapping_rate = reads.empty() ? 0.0 : ((double)mappingCount / reads.size() * 100.0);

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n걸린 시간              : " << elapsed_sec         << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory              << " MB\n";
    cout << "총 원본 염기서열 길이(N): " << N                   << "\n";
    cout << "일치하지 않는 글자 수   : " << mismatched          << " 개 (미복구 빈칸 포함)\n";
    cout << "재구축 일치율(정확도)   : " << reconstruction_rate << "%\n";

    // 재구축된 염기서열 파일 저장
    ofstream fout("reconstructed_greedy_seq.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_greedy\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_greedy_seq.txt\n";
    } else {
        cerr << "Error: reconstructed_greedy_seq.txt 저장 실패\n";
    }

    // =========================================================================
    // [CSV 저장 파트] 결과 수집 후 아래 블록 전체 삭제 예정
    // =========================================================================
    bool is_new_file = false;
    ifstream check_file("results_redblackTree_Greedy.csv");
    if (!check_file.is_open()) {
        is_new_file = true;
    } else {
        check_file.close();
    }

    ofstream csv("results_redblackTree_Greedy.csv", ios::app);
    if (csv.is_open()) {
        if (is_new_file) {
            csv << "algorithm,dataset,snp_rate,total_sec,memory_mb,mapped_pct,reconstruct_pct\n";
        }
        csv << fixed << setprecision(2)
            << current_algorithm << "," << current_dataset << "," << current_snp << ","
            << elapsed_sec << "," << memory << "," << mapping_rate << ","
            << reconstruction_rate << "\n";
        csv.close();
        cout << "\n결과가 'results_redblackTree_Greedy.csv' 에 기록되었습니다!\n";
    } else {
        cout << "\nCSV 파일 저장 실패\n";
    }
    // =========================================================================

    return 0;
}