#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

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

    string result(N, 'N');
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

// ── 정확도 측정 ──
double calcAccuracy(const string& original, const string& assembled) {
    int match = 0;
    int n = min(original.size(), assembled.size());

    for (int i = 0; i < n; i++) {
        if (assembled[i] != 'N' && original[i] == assembled[i]) {
            match++;
        }
    }

    return (double)match / original.size() * 100.0;
}

// ── main ──
/*int main() {
    auto startTime = chrono::high_resolution_clock::now();

    // 1. 파일 읽기
    // [1] 인공 서열
    string ref         = load_reference("reference_synthetic.txt");
    vector<Read> reads = load_reads("reads_synthetic.txt");
    string original    = load_original("original_synthetic_1M.txt");

    // [2] 빵효모 (반복서열 영향 확인)
    //string ref         = load_reference("reference_yeast.txt");
    //vector<Read> reads = load_reads("reads_yeast.txt");
    //string original    = load_original("original_yeast_1M.txt");

    // [3] Baseline (에러 없는 깨끗한 read)
    //string ref         = load_reference(reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_baseline.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    // [4] InDel (삽입/결실 폭격 - DP 가 빛나는 환경)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_indel.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    // [5] End-Heavy (read 후반부에 에러 집중)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_end_heavy.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    if (ref.empty() || reads.empty()) {
        cout << "데이터 로딩 실패" << endl;
        return 1;
    }

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
    string assembled = consensus(ref, reads, positions);

    // 5. 정확도 측정
    double accuracy = calcAccuracy(original, assembled);

    // 6. 시간 측정
    auto endTime = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(endTime - startTime).count();

    // 7. 메모리 사용량 측정
    double sa_memory_mb = (double)(sa.size() * sizeof(int)) / (1024.0 * 1024.0);

    // 8. 결과 출력
    cout.setf(ios::fixed);

    cout << "\n=== 결과 ===" << endl;
    cout << "총 수행 시간: " << elapsed << "초" << endl;
    cout << "매핑된 read: " << mapped << " / " << reads.size() << endl;

    // reads.size()가 0일 때의 예외 처리 추가
    double mapping_rate = reads.empty() ? 0.0 : ((double)mapped / reads.size() * 100.0);
    cout << "매핑률: " << mapping_rate << "%" << endl;

    cout << "복원 정확도: " << accuracy << "%" << endl;
    cout << "SA 메모리 사용량: " << sa_memory_mb << " MB" << endl;

    // 8. 복원된 염기서열 파일로 저장 (FASTA 포맷)
    ofstream fout("reconstructed_sequence_SA.fasta");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_SA_mapped\n";

        // 생물정보학 표준 FASTA 형식에 맞게 80글자마다 줄바꿈을 해줍니다.
        for (size_t i = 0; i < assembled.size(); i += 80) {
            fout << assembled.substr(i, 80) << "\n";
        }
        fout.close();
    }

    return 0;
}*/
int main() {
    // =========================================================================
    // [사용자 설정 파트]
    // 코드를 실행하기 전에 현재 돌리는 데이터셋의 정보를 여기에 적어주세요!
    // 이 정보가 그대로 CSV 파일에 기록됩니다.
    // =========================================================================
    string current_algorithm = "sa_minimizer";
    string current_snp       = "0.1%";       // 예: 0.1%, 1.0%, 2.0%, 5.0%
    string current_dataset   = "yeast";  // 예: synthetic, baseline, indel, end_heavy
    // =========================================================================

    // 1. 파일 읽기 (수빈님이 하시던 대로 주석 풀고 닫기 하시면 됩니다)
    // [1] 인공 서열
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_synthetic.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    // [2] 빵효모 (반복서열 영향 확인)
    //string ref         = load_reference("reference_yeast.txt");
    //vector<Read> reads = load_reads("reads_yeast.txt");
    //string original    = load_original("original_yeast_1M.txt");

    // [3] Baseline (에러 없는 깨끗한 read)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_baseline.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    // [4] InDel (삽입/결실 폭격 - DP 가 빛나는 환경)
    string ref         = load_reference("reference_synthetic.txt");
    vector<Read> reads = load_reads("reads_indel.txt");
    string original    = load_original("original_synthetic_1M.txt");

    // [5] End-Heavy (read 후반부에 에러 집중)
    //string ref         = load_reference("reference_synthetic.txt");
    //vector<Read> reads = load_reads("reads_end_heavy.txt");
    //string original    = load_original("original_synthetic_1M.txt");

    if (ref.empty() || reads.empty()) {
        cout << "데이터 로딩 실패" << endl;
        return 1;
    }

    // 2. Suffix Array 구축 (시간 측정 분리)
    auto idx_start = chrono::high_resolution_clock::now();
    vector<int> sa = buildSuffixArray(ref);
    auto idx_end = chrono::high_resolution_clock::now();
    double index_sec = chrono::duration<double>(idx_end - idx_start).count();

    // 3. Read 매핑 (Search 시간 측정 분리)
    int k = 15;
    int max_mismatch = 2;

    vector<int> positions(reads.size(), -1);
    int mapped = 0;
    int correct_position = 0; // 위치 정확도 계산용 변수

    auto srch_start = chrono::high_resolution_clock::now();
    for (int i = 0; i < (int)reads.size(); i++) {
        positions[i] = mapRead(sa, ref, reads[i].sequence, k, max_mismatch);
        if (positions[i] != -1) {
            mapped++;
            // 매핑된 위치가 정답 위치(start_pos)와 완벽히 일치하는지 체크
            if (positions[i] == reads[i].start_pos) {
                correct_position++;
            }
        }
    }
    auto srch_end = chrono::high_resolution_clock::now();
    double search_sec = chrono::duration<double>(srch_end - srch_start).count();

    double total_sec = index_sec + search_sec;

    // 4. 서열 복원
    string assembled = consensus(ref, reads, positions);

    // 5. 정확도 측정
    double accuracy = calcAccuracy(original, assembled);

    // 6. 메모리 사용량 및 통계 계산
    double sa_memory_mb = (double)(sa.size() * sizeof(int)) / (1024.0 * 1024.0);
    double mapping_rate = reads.empty() ? 0.0 : ((double)mapped / reads.size() * 100.0);
    double position_accuracy_pct = reads.empty() ? 0.0 : ((double)correct_position / reads.size() * 100.0);

    // 7. 결과 출력 (콘솔)
    cout.setf(ios::fixed);
    cout << "\n=== 결과 ===" << endl;
    cout << "총 수행 시간: " << total_sec << "초" << endl;
    cout << "매핑된 read: " << mapped << " / " << reads.size() << endl;
    cout << "매핑률: " << mapping_rate << "%" << endl;
    cout << "복원 정확도: " << accuracy << "%" << endl;
    cout << "SA 메모리 사용량: " << sa_memory_mb << " MB" << endl;

    // 8. 복원된 염기서열 파일로 저장
    ofstream fout("reconstructed_sequence_SA.fasta");
    if (fout.is_open()) {
        fout << ">reconstructed_sequence_SA_mapped\n";
        for (size_t i = 0; i < assembled.size(); i += 80) {
            fout << assembled.substr(i, 80) << "\n";
        }
        fout.close();
    }

    // 9. CSV 파일에 누적 저장 (Append Mode)
    bool is_new_file = false;
    ifstream check_file("results_sa_minimizer.csv");
    if (!check_file.is_open()) {
        is_new_file = true; // 파일이 없으면 새로 생성하며 헤더를 적기 위해 체크
    } else {
        check_file.close();
    }

    // ios::app 모드로 열면 기존 파일 내용을 지우지 않고 맨 아래에 이어 씁니다.
    ofstream csv("results_sa_minimizer.csv", ios::app);
    if (csv.is_open()) {
        // 파일이 처음 만들어질 때만 첫 줄(헤더) 작성
        if (is_new_file) {
            csv << "algorithm,snp_rate,dataset,index_sec,search_sec,total_sec,memory_mb,mapped_pct,position_accuracy_pct,reconstruct_pct\n";
        }
        // 현재 실행된 결과 데이터를 한 줄 추가
        csv << fixed << setprecision(2)
            << current_algorithm << "," << current_snp << "," << current_dataset << ","
            << index_sec << "," << search_sec << "," << total_sec << ","
            << sa_memory_mb << "," << mapping_rate << ","
            << position_accuracy_pct << "," << accuracy << "\n";
        csv.close();
        cout << "\n✅ 'results_sa_minimizer.csv' 파일에 결과가 성공적으로 기록되었습니다!" << endl;
    } else {
        cout << "\n❌ CSV 파일 저장 실패" << endl;
    }

    return 0;
}