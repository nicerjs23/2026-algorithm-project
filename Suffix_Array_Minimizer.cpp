#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sys/time.h>
#include <sys/resource.h>
#include <iomanip>

using namespace std;

// ── 구조체 ──
struct Read {
    int id;
    int start_pos;
    string sequence;
};

struct MinimizerInfo {
    string sequence;
    int pos_in_read;
};

// ── 파일 로드 ──
string load_reference(string filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "에러: " << filename << " 을 찾을 수 없습니다." << endl;
        exit(1);
    }
    string ref;
    fin >> ref;
    return ref;
}

vector<Read> load_reads(string filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "에러: " << filename << " 을 찾을 수 없습니다." << endl;
        exit(1);
    }
    vector<Read> reads;
    string header;
    getline(fin, header); // 헤더 스킵
    Read r;
    while (fin >> r.id >> r.start_pos >> r.sequence) {
        reads.push_back(r);
    }
    return reads;
}

string load_original(string filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cout << "에러: " << filename << " 을 찾을 수 없습니다." << endl;
        exit(1);
    }
    string orig;
    fin >> orig;
    return orig;
}

// ── Suffix Array 구축 ──
bool compare_suffixes(int a, int b, const string& text) {
    int n = text.length();
    while (a < n && b < n) {
        if (text[a] != text[b]) return text[a] < text[b];
        a++; b++;
    }
    return a == n;
}

int partition(vector<int>& sa, int left, int right, const string& text) {
    int mid = left + (right - 1 - left) / 2;
    swap(sa[left], sa[mid]);

    int pivot_idx = left;
    int pivot = sa[pivot_idx];

    do {
        do { left++; } while (left < right && compare_suffixes(sa[left], pivot, text));
        do { right--; } while (right > pivot_idx && compare_suffixes(pivot, sa[right], text));
        if (left < right) swap(sa[left], sa[right]);
        else break;
    } while (true);

    sa[pivot_idx] = sa[right];
    sa[right] = pivot;
    return right;
}

void iterativeQuickSort(vector<int>& sa, int left, int right, const string& text) {
    vector<int> stack;
    stack.push_back(left);
    stack.push_back(right);

    while (!stack.empty()) {
        int r = stack.back(); stack.pop_back();
        int l = stack.back(); stack.pop_back();

        if (l < r) {
            int k = partition(sa, l, r + 1, text);
            stack.push_back(l);   stack.push_back(k - 1);
            stack.push_back(k + 1); stack.push_back(r);
        }
    }
}

vector<int> buildSuffixArray(const string& text) {
    int n = text.length();
    vector<int> sa(n);
    for (int i = 0; i < n; i++) sa[i] = i;
    iterativeQuickSort(sa, 0, n - 1, text);
    return sa;
}

// ── Minimizer 추출 ──
MinimizerInfo extractMinimizer(const string& read_seq, int k) {
    MinimizerInfo best;
    best.sequence = read_seq.substr(0, k);
    best.pos_in_read = 0;

    for (int i = 1; i <= (int)read_seq.length() - k; i++) {
        string kmer = read_seq.substr(i, k);
        if (kmer < best.sequence) {
            best.sequence = kmer;
            best.pos_in_read = i;
        }
    }
    return best;
}

// ── SA 이진탐색 ──
int searchSA(const vector<int>& sa, const string& ref, const string& target) {
    int n = sa.size();
    int k = target.length();
    int left = 0, right = n - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        int pos = sa[mid];
        string kmer = ref.substr(pos, min(k, n - pos));

        if (kmer == target) return pos;
        if (target < kmer) right = mid - 1;
        else left = mid + 1;
    }
    return -1;
}

// ── mismatch 계산 ──
int countMismatch(const string& ref, const string& read_seq, int start_pos) {
    if (start_pos < 0 || start_pos + (int)read_seq.length() > (int)ref.length())
        return 999;

    int count = 0;
    for (int i = 0; i < (int)read_seq.length(); i++) {
        if (ref[start_pos + i] != read_seq[i]) count++;
    }
    return count;
}

int mapRead(const vector<int>& sa, const string& ref,
            const string& read_seq, int k, int max_mismatch) {

    // Minimizer 추출
    MinimizerInfo min_info = extractMinimizer(read_seq, k);

    // SA에서 위치 탐색
    int found_pos = searchSA(sa, ref, min_info.sequence);
    if (found_pos == -1) return -1;

    // 실제 read 시작 위치 계산
    int real_start = found_pos - min_info.pos_in_read;

    // mismatch 검증
    int mm = countMismatch(ref, read_seq, real_start);
    if (mm <= max_mismatch) return real_start;

    return -1;
}

// ── 투표로 서열 복원 ──
string consensus(const string& ref,
                 const vector<Read>& reads,
                 const vector<int>& positions) {
    int N = ref.size();
    vector<array<int,4>> votes(N, {0,0,0,0});
    string bases = "ACGT";

    for (int i = 0; i < (int)reads.size(); i++) {
        int pos = positions[i];
        if (pos < 0) continue;

        for (int j = 0; j < (int)reads[i].sequence.size(); j++) {
            if (pos + j >= N) break;
            char base = reads[i].sequence[j];
            if      (base == 'A') votes[pos+j][0]++;
            else if (base == 'C') votes[pos+j][1]++;
            else if (base == 'G') votes[pos+j][2]++;
            else if (base == 'T') votes[pos+j][3]++;
        }
    }

    string result(N, '-');
    for (int i = 0; i < N; i++) {
        int maxVote = 0, maxIdx = 0;
        for (int j = 0; j < 4; j++) {
            if (votes[i][j] > maxVote) {
                maxVote = votes[i][j];
                maxIdx = j;
            }
        }
        if (maxVote > 0) result[i] = bases[maxIdx];
    }
    return result;
}

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

// ── main ──
int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);  // Windows 콘솔 UTF-8
#endif

    // =========================================================================
    // [사용자 설정 파트]
    // =========================================================================
    string current_algorithm = "sa_minimizer";
    string current_snp       = "1.0%";       // 예: 0.1%, 1.0%, 2.0%, 5.0%
    string current_dataset   = "baseline";   // 예: yeast, baseline, indel, end_heavy
    // =========================================================================

    cout << "데이터 로딩 중...\n";

    // 1. 파일 읽기
    // [1] Baseline (에러 없는 깨끗한 read) - 인공 서열
    string ref         = load_reference("reference_synthetic.txt");
    vector<Read> reads = load_reads("reads_baseline.txt");
    string original    = load_original("original_synthetic_1M.txt");

    // [2] 빵효모 (반복서열 영향 확인)
    //string ref         = load_reference("reference_yeast.txt");
    //vector<Read> reads = load_reads("reads_yeast.txt");
    //string original    = load_original("original_yeast_1M.txt");

    // [3] InDel (삽입/결실 폭격 - DP 가 빛나는 환경)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_indel.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    // [4] End-Heavy (read 후반부에 에러 집중)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_end_heavy.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    if (ref.empty() || reads.empty()) {
        cout << "데이터 로딩 실패" << endl;
        return 1;
    }
    cout << "SA_Minimizer mapping...\n";

    // ── 알고리즘 순수 성능 시간 측정 시작 ──
    clock_t start_time = clock();

    // 2. Suffix Array 구축
    vector<int> sa = buildSuffixArray(ref);

    // 3. Read 매핑
    int k = 15;
    int max_mismatch = 2;
    vector<int> positions(reads.size(), -1);
    int mapped = 0;

    for (int i = 0; i < (int)reads.size(); i++) {
        positions[i] = mapRead(sa, ref, reads[i].sequence, k, max_mismatch);
        if (positions[i] != -1) mapped++;
    }

    // 4. 서열 복원
    string reconstructed_seq = consensus(ref, reads, positions);

    // 5. 정확도 측정 (팀원 코드와 100% 동일한 로직 적용)
    int N = (int)ref.length();
    int mismatched = 0;
    for (int i = 0; i < N; ++i) {
        if (reconstructed_seq[i] != original[i]) {
            ++mismatched;
        }
    }

    clock_t end_time = clock();

    // 계산
    double elapsed_sec = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    double sa_memory_mb = check_memory();
    double correct_rate = 100.0 * (N - mismatched) / N;
    double mapping_rate = reads.empty() ? 0.0 : ((double)mapped / reads.size() * 100.0);

    // 8. 결과 출력
    cout.setf(ios::fixed);
    cout.precision(2);

    cout << "\n=== 결과 ===" << endl;
    cout << "총 수행 시간: " << elapsed_sec << " 초\n";
    cout << "매핑된 read: " << mapped << " / " << reads.size() << endl;
    cout << "매핑률: " << mapping_rate << "%\n";
    cout << "복원 정확도: " << correct_rate << "%\n";
    cout << "SA 메모리 사용량: " << sa_memory_mb << " MB\n";

    // 9. 복원된 염기서열 파일로 저장
    ofstream fout("reconstructed_sequence_SA.fasta");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_SA_mapped\n";
        for (int i = 0; i < N; i += 60) {
            fout << reconstructed_seq.substr(i, 60) << "\n";
        }
        fout.close();
        cout << "재구축 서열 저장 완료  : reconstructed_sequence_SA.fasta\n";
    }

    // 10. CSV 파일에 누적 저장
    bool is_new_file = false;
    ifstream check_file("results_sa_minimizer.csv");
    if (!check_file.is_open()) {
        is_new_file = true;
    } else {
        check_file.close();
    }

    ofstream csv("results_sa_minimizer.csv", ios::app);
    if (csv.is_open()) {
        if (is_new_file) {
            csv << "algorithm,dataset,snp_rate,total_sec,memory_mb,mapped_pct,reconstruct_pct\n";
        }
        csv << fixed << setprecision(2)
            << current_algorithm << "," << current_dataset << "," << current_snp << ","
            << elapsed_sec << "," << sa_memory_mb << "," << mapping_rate << ","
            << correct_rate << "\n";
        csv.close();
        cout << "\n✅ 핵심 데이터가 'results_sa_minimizer.csv' 파일에 성공적으로 기록되었습니다!\n";
    } else {
        cout << "\n❌ CSV 파일 저장 실패\n";
    }

    return 0;
}
