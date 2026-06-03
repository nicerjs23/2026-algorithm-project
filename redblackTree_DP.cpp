#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <map>
#include <algorithm>
#include <windows.h>
#include <psapi.h>

using namespace std;

// ============================================================
//  상수 정의 및 환경 설정
// ============================================================
const int K_MER = 10;
const int READ_LEN = 30;

// DP 점수 체계
const int SCORE_MATCH = 2;
const int SCORE_MISMATCH = -1;
const int SCORE_GAP = -2;

// 오타를 최대 2개까지만 허용하는 논리적 최소 임계값 자동 계산
const int SCORE_THRESHOLD = (READ_LEN * SCORE_MATCH) + 2 * (SCORE_MISMATCH - SCORE_MATCH);

// ============================================================
//  구조체
// ============================================================
struct Read {
    int start_pos;
    string seq;
};

// 메모리 측정 (실제 프로세스 메모리 - psapi)
double check_memory() {
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return (double)pmc.WorkingSetSize / (1024.0 * 1024.0);
}

// 레퍼런스 및 원본 파일 읽기
string readReference(const string& filename) {
    ifstream input(filename);
    string content;
    if (input.is_open()) getline(input, content);
    return content;
}

// Reads 파일 읽기 (탭 구분자 파싱)
vector<Read> readReads(const string& filename) {
    ifstream input(filename);
    vector<Read> reads;
    if (!input.is_open()) return reads;
    reads.reserve(100000);
    string line;
    getline(input, line); // 헤더 스킵
    while (getline(input, line)) {
        size_t tab1 = line.find('\t');
        size_t tab2 = line.find('\t', tab1 + 1);
        Read r;
        r.start_pos = stoi(line.substr(tab1 + 1, tab2 - tab1 - 1));
        r.seq = line.substr(tab2 + 1);
        reads.push_back(r);
    }
    return reads;
}

// ============================================================
// [1단계] Red-Black Tree 기반 Seed 색인 구축
// ============================================================
map<string, vector<int>> buildSeedIndex(const string& reference, int k) {
    map<string, vector<int>> index;
    int N = (int)reference.length();
    for (int i = 0; i <= N - k; ++i) {
        index[reference.substr(i, k)].push_back(i);
    }
    return index;
}

// ============================================================
// [2단계] 메모리 초고도화 최적화 (1D 배열 기반 DP 계산)
// ============================================================
int calculateDPScore(const string& read_seq, const string& ref_window) {
    int n = read_seq.length();
    int m = ref_window.length();

    vector<int> prev_row(m + 1, 0);
    vector<int> curr_row(m + 1, 0);

    for (int j = 1; j <= m; j++) prev_row[j] = j * SCORE_GAP;

    for (int i = 1; i <= n; i++) {
        curr_row[0] = i * SCORE_GAP;
        for (int j = 1; j <= m; j++) {
            int scoreDiag = prev_row[j - 1] + (read_seq[i - 1] == ref_window[j - 1] ? SCORE_MATCH : SCORE_MISMATCH);
            int scoreUp = prev_row[j] + SCORE_GAP;
            int scoreLeft = curr_row[j - 1] + SCORE_GAP;
            curr_row[j] = max({scoreDiag, scoreUp, scoreLeft});
        }
        prev_row = curr_row;
    }

    int best = prev_row[0];
    for (int j = 1; j <= m; j++) {
        best = max(best, prev_row[j]);
    }
    return best;
}

// ============================================================
//  메인 제어 흐름
// ============================================================
int main() {
    SetConsoleOutputCP(65001);

    cout << "데이터 로딩 중...\n";
    
    // 인공서열
    string reference_genome = readReference("reference_synthetic.txt");
    vector<Read> reads      = readReads("reads_baseline.txt");
    // vector<Read> reads = readReads("reads_indel.txt");
    // vector<Read> reads = readReads("reads_end_heavy.txt");
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
    string reconstructed_seq(N, '-');

    cout << "Index 테이블 생성 중...\n";
    map<string, vector<int>> seed_index = buildSeedIndex(reference_genome, K_MER);

    cout << "Seed-and-Extend 매핑 시작...\n";
    clock_t start_time = clock();

    int mapped_count = 0;
    int correct_mapping_count = 0;

    for (size_t i = 0; i < reads.size(); ++i) {
        int best_pos = -1;
        int best_score = SCORE_THRESHOLD - 1;

        for (int offset = 0; offset <= 20; offset += 10) {
            string seed = reads[i].seq.substr(offset, K_MER);
            auto it = seed_index.find(seed);

            if (it != seed_index.end()) {
                for (int pos : it->second) {
                    int ref_start = pos - offset;
                    if (ref_start < 0 || ref_start >= N) continue;

                    int window_len = min(READ_LEN + 5, N - ref_start);
                    string ref_window = reference_genome.substr(ref_start, window_len);

                    int score = calculateDPScore(reads[i].seq, ref_window);

                    if (score > best_score) {
                        best_score = score;
                        best_pos = ref_start;
                    }
                }
            }
        }

        if (best_pos != -1) {
            mapped_count++;
            if (best_pos == reads[i].start_pos) {
                correct_mapping_count++;
            }
            for (int j = 0; j < READ_LEN; ++j) {
                if (best_pos + j < N) {
                    reconstructed_seq[best_pos + j] = reads[i].seq[j];
                }
            }
        }
    }

    // 최종 검증
    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original_seq[i]) ++mismatched;
    }

    clock_t end_time = clock();
    double elapsed_sec       = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double memory            = check_memory();  // 실제 프로세스 메모리
    double reconstruction_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n=========================================\n";
    cout << "걸린 시간              : " << elapsed_sec         << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory              << " MB\n";
    cout << "총 원본 염기서열 길이(N): " << N                   << "\n";
    cout << "일치하지 않는 글자 수   : " << mismatched          << " 개 (미복구 빈칸 포함)\n";
    cout << "재구축 일치율(정확도)   : " << reconstruction_rate << "%\n";
    cout << "=========================================\n";

    // 재구축된 염기서열 파일 저장 (FASTA 형식, 60bp마다 줄바꿈)
    ofstream fout("reconstructed_seq.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed_seq.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_seq.txt\n";
    } else {
        cerr << "Error: reconstructed_seq.txt 저장 실패\n";
    }

    return 0;
}