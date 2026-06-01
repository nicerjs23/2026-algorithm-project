#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <windows.h>

using namespace std;

// ============================================================
//  구조체 정의
// ============================================================

struct Read {
    int start_pos;
    string seq;
};

struct Result {
    int total = 0;
    int mapped = 0;
    int correct = 0;
};

// ============================================================
//  메모리 측정
// ============================================================
double check_memory(const string& ref, const vector<Read>& reads) {
    size_t mem_bytes = 0;

    mem_bytes += ref.capacity();
    mem_bytes += reads.capacity() * sizeof(Read);
    for(size_t i = 0; i < reads.size(); ++i) {
        mem_bytes += reads[i].seq.capacity();
    }

    return (double)mem_bytes / (1024.0 * 1024.0);
}

// ============================================================
//  파일 로딩
// ============================================================

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

    getline(input, line); // 1. 헤더 스킵

    // 2. 파일 끝까지 한 줄씩 읽기
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
// Brute-force Search
// ============================================================

int bruteForceSearch(const string& text, const string& pattern, int mismatch = 2) {
    const int N = (int)text.size();
    const int L = (int)pattern.size();

    int pos = -1;
    int cut  = mismatch + 1;

    for (int i = 0; i <= N - L; ++i) {
        int mm = 0;
        for (int j = 0; j < L; ++j) {
            if (text[i + j] != pattern[j]) {
                ++mm;
                if (mm > mismatch) break;
            }
        }
        if (mm <= mismatch && mm < cut) {
            cut = mm;
            pos = i;
        }

        if (cut == 0) break; // 완벽 매칭이면 즉시 종료
    }
    return pos; // 위치 하나만 깔끔하게 반환!
}

// ============================================================
//  메인
// ============================================================

int main() {
    SetConsoleOutputCP(65001);
    cout << "데이터 로딩 중...\n";

    clock_t start_time = clock();

    string reference_genome = readReference("reference_genome.txt");
    vector<Read> reads      = readReads("reads.txt");

    string original_seq = readReference("original_1M.txt");

    if (reference_genome.empty() || reads.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    cout << "brute force mapping...\n";

    int N = reference_genome.length();
    string reconstructed_seq(N, '-');

    Result result;
    result.total = (int)reads.size();

    int progress_step = result.total / 10;
    for (int i = 0; i < result.total; ++i) {

        int pos = bruteForceSearch(reference_genome, reads[i].seq);

        if (pos != -1) {
            ++result.mapped;

            // Reconstruct
            int L = reads[i].seq.length();
            for (int j = 0; j < L; ++j) {

                if (pos + j < N) {
                    reconstructed_seq[pos + j] = reads[i].seq[j];
                }
            }
        }

        if (progress_step > 0 && (i + 1) % progress_step == 0) {
            cout << "  Progress: " << (i + 1) << " / " << result.total
                 << " (" << (100 * (i + 1) / result.total) << "%)\n";
        }
    }

    int mismatched = 0; // 일치하지 않는 a,c,g,t 개수
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original_seq[i]) {
            ++mismatched;
        }
    }

    clock_t end_time = clock();
    double elapsed_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double memory = check_memory(reference_genome, reads);
    double correct_rate = 100.0 * (N - mismatched) / N; // 글자 일치율

    cout.setf(ios::fixed);
    cout.precision(2);

    // 최종 결과 출력
    cout << "\n================================\n";
    cout << "걸린 시간              : " << elapsed_sec << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory << " MB\n";
    cout << "총 원본 염기서열 길이         : " << N << "\n";
    cout << "일치하지 않는 글자 수           : " << mismatched << "\n";
    cout << "재구축 일치율                 : " << correct_rate << "%\n";
    cout << "================================\n";

    return 0;
}