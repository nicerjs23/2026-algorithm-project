#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <map>
#include <algorithm>

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

// DP 점수 체계
const int dpMatch = 2;
const int dpMismatch = -1;
const int dpIndel = -2;

// mismatch 임계값 계산
const int isMatch = (readLen * dpMatch) + 2 * (dpMismatch - dpMatch);

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
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);  // macOS: bytes
  #else
    return (double)ru.ru_maxrss / 1024.0;             // Linux: KB
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

// DP 계산
int calculateDPScore(const string& readSeq, const string& refCandidate) {
    int n = readSeq.length();
    int m = refCandidate.length();

    vector<int> prevRow(m + 1, 0);
    vector<int> curRow(m + 1, 0);

    for (int j = 1; j <= m; j++) {
        prevRow[j] = j * dpIndel;
    }

    for (int i = 1; i <= n; i++) {
        curRow[0] = i * dpIndel;
        for (int j = 1; j <= m; j++) {
            int scoreDiag = prevRow[j - 1] + (readSeq[i - 1] == refCandidate[j - 1] ? dpMatch : dpMismatch); // snp인 경우만 고려
            int scoreUp = prevRow[j] + dpIndel; // 결실 발생
            int scoreLeft = curRow[j - 1] + dpIndel; // 삽입 발생
            curRow[j] = max({scoreDiag, scoreUp, scoreLeft}); // 세 가지 중 점수가 가장 높은 경우 선택
        }
        prevRow = curRow;
    }

    int best = prevRow[0];
    for (int j = 1; j <= m; j++) {
        best = max(best, prevRow[j]);
    }
    return best;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8 (한글 깨짐 방지)
#endif

    cout << "데이터 로딩 중...\n";

    // 인공서열
    string reference_genome = readReference("reference_synthetic.txt");
    //vector<Read> reads      = readReads("reads_baseline.txt");
    //vector<Read> reads = readReads("reads_indel.txt");
    vector<Read> reads = readReads("reads_end_heavy.txt");
    string original_seq     = readReference("original_synthetic_1M.txt");

    // 효모
    //string reference_genome = readReference("reference_yeast.txt");
    //vector<Read> reads      = readReads("reads_yeast.txt");
    //string original_seq     = readReference("original_yeast_1M.txt");

    if (reference_genome.empty() || reads.empty() || original_seq.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    int N = reference_genome.length();
    string reconstructed(N, '-');

    cout << "Red-Black Tree Seed-and-Extend...\n";
    map<string, vector<int>> seedIdx = consSeedIdx(reference_genome, k);

    cout << "Seed-and-Extend 매핑...\n";
    clock_t sTime = clock();

    int mappingCount = 0;
    int exactCount = 0;

    for (size_t i = 0; i < reads.size(); ++i) {
        int bestPos = -1;
        int bestScore = isMatch;

        for (int offset = 0; offset <= 20; offset += 10) {
            string seed = reads[i].seq.substr(offset, k);
            auto it = seedIdx.find(seed); // RB-tree에서 mismatch가 없는 구간 찾기 (find)

            if (it != seedIdx.end()) {
                for (int pos : it->second) {
                    int refStart = pos - offset;
                    if (refStart < 0 || refStart >= N) continue;

                    int candidateLen = min(readLen + 2, N - refStart);
                    string refCandidate = reference_genome.substr(refStart, candidateLen);

                    int score = calculateDPScore(reads[i].seq, refCandidate);

                    if (score >= bestScore) {
                        bestScore = score;
                        bestPos = refStart;
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
    double memory = checkMemory();  // 실제 프로세스 메모리
    double reconstruction_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n걸린 시간              : " << elapsed_sec         << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory              << " MB\n";
    cout << "총 원본 염기서열 길이(N): " << N                   << "\n";
    cout << "일치하지 않는 글자 수   : " << mismatched          << " 개 (미복구 빈칸 포함)\n";
    cout << "재구축 일치율(정확도)   : " << reconstruction_rate << "%\n";

    // 재구축된 염기서열 파일 저장
    ofstream fout("reconstructed_seq.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_seq.txt\n";
    } else {
        cerr << "Error: reconstructed_seq.txt 저장 실패\n";
    }

    return 0;
}