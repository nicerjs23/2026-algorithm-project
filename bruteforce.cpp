#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>

// ---- 플랫폼별 메모리 측정 / 콘솔 인코딩 분기 ----
#ifdef _WIN32
  #include <windows.h>
  #include <psapi.h>
#else
  #include <sys/resource.h>
#endif

using namespace std;

struct Read {
    int start_pos;
    string seq;
};

struct Result {
    int total  = 0;
    int mapped = 0;
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
    getline(input, line); // 헤더 스킵
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

int bruteForceSearch(const string& text, const string& pattern, int mismatch = 2) {
    const int N = (int)text.size();
    const int L = (int)pattern.size();

    int pos = -1;
    int cut = mismatch + 1;

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
        if (cut == 0) break;
    }
    return pos;
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8
#endif
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

    if (reference_genome.empty() || reads.empty()) {
        cerr << "Error: 데이터 로딩 실패\n";
        return 1;
    }

    cout << "brute force mapping...\n";
    clock_t start_time = clock();

    int N = (int)reference_genome.length();
    string reconstructed_seq(N, '-');

    Result result;
    result.total = (int)reads.size();

    int progress_step = result.total / 10;
    for (int i = 0; i < result.total; ++i) {
        int pos = bruteForceSearch(reference_genome, reads[i].seq);

        if (pos != -1) {
            ++result.mapped;

            int L = (int)reads[i].seq.length();
            for (int j = 0; j < L; ++j) {
                if (pos + j < N)
                    reconstructed_seq[pos + j] = reads[i].seq[j];
            }
        }

        if (progress_step > 0 && (i + 1) % progress_step == 0) {
            cout << "  Progress: " << (i + 1) << " / " << result.total
                 << " (" << (100 * (i + 1) / result.total) << "%)\n";
        }
    }

    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original_seq[i])
            ++mismatched;
    }

    clock_t end_time = clock();
    double elapsed_sec  = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double memory       = check_memory();
    double correct_rate = 100.0 * (N - mismatched) / N;

    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n걸린 시간              : " << elapsed_sec  << " 초\n";
    cout << "사용 중인 메모리 크기   : " << memory       << " MB\n";
    cout << "총 원본 염기서열 길이   : " << N            << "\n";
    cout << "일치하지 않는 글자 수   : " << mismatched   << "\n";
    cout << "재구축 일치율           : " << correct_rate << "%\n";

    // 재구축된 염기서열 파일 저장
    ofstream fout("reconstructedbf_seq.txt");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed_seq.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructedbf_seq.txt\n";
    } else {
        cerr << "Error: reconstructedbf_seq.txt 저장 실패\n";
    }

    return 0;
}