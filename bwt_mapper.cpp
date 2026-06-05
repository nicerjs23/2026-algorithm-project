// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할2 (김세훈)
//  BWT + FM-index (Suffix Array + Backward Search)
//  16가지 실험 자동 실행 → results_bwt.csv 저장
//
//  [실험 조합]
//   Synthetic: baseline / indel / end_heavy × 0.1% / 1.0% / 2.0% / 5.0% = 12
//   Yeast    : yeast × 0.1% / 1.0% / 2.0% / 5.0% = 4
//   합계 = 16
//
//  [최적화]
//   같은 reference는 인덱스를 한 번만 구축 (8번 구축으로 16번 실험)
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <ctime>

#ifdef _WIN32
  #include <windows.h>
  #include <psapi.h>
#else
  #include <sys/resource.h>
#endif

using namespace std;

static const int SEED_LEN     = 10;
static const int MAX_MISMATCH = 2;
static const int MAX_CAND     = 200;
static const int ALPHA        = 5;

struct Read { int start_pos; string seq; };

double check_memory() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return (double)pmc.PeakWorkingSetSize / (1024.0 * 1024.0);
#else
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
  #ifdef __APPLE__
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
  #else
    return (double)ru.ru_maxrss / 1024.0;
  #endif
#endif
}

static int cidx(char c) {
    switch (c) {
        case '$': return 0; case 'A': return 1;
        case 'C': return 2; case 'G': return 3;
        case 'T': return 4; default:  return -1;
    }
}

vector<int> buildSA(const string& s) {
    int n = (int)s.size();
    vector<int> sa(n), rank_(n), tmp(n);
    iota(sa.begin(), sa.end(), 0);
    for (int i = 0; i < n; i++) rank_[i] = (unsigned char)s[i];
    for (int gap = 1; gap < n; gap <<= 1) {
        sort(sa.begin(), sa.end(), [&](int a, int b) {
            if (rank_[a] != rank_[b]) return rank_[a] < rank_[b];
            int ra = (a+gap<n)?rank_[a+gap]:-1;
            int rb = (b+gap<n)?rank_[b+gap]:-1;
            return ra < rb;
        });
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) {
            tmp[sa[i]] = tmp[sa[i-1]];
            int ra1=rank_[sa[i]], rb1=rank_[sa[i-1]];
            int ra2=(sa[i]+gap<n)?rank_[sa[i]+gap]:-1;
            int rb2=(sa[i-1]+gap<n)?rank_[sa[i-1]+gap]:-1;
            if (ra1!=rb1||ra2!=rb2) tmp[sa[i]]++;
        }
        rank_ = tmp;
        if (rank_[sa[n-1]] == n-1) break;
    }
    return sa;
}

string buildBWT(const string& s, const vector<int>& sa) {
    int n = (int)s.size();
    string bwt(n, ' ');
    for (int i = 0; i < n; i++)
        bwt[i] = (sa[i]==0) ? '$' : s[sa[i]-1];
    return bwt;
}

struct FMIndex {
    vector<int> C;
    vector<array<int,ALPHA>> Occ;
    explicit FMIndex(const string& bwt) {
        int n = (int)bwt.size();
        Occ.assign(n+1, {});
        for (int i = 0; i < n; i++) {
            Occ[i+1] = Occ[i];
            int c = cidx(bwt[i]);
            if (c >= 0) Occ[i+1][c]++;
        }
        C.resize(ALPHA, 0);
        for (int c = 1; c < ALPHA; c++)
            C[c] = C[c-1] + Occ[n][c-1];
    }
};

pair<int,int> backwardSearch(const string& pat, const FMIndex& fm) {
    int n = (int)fm.Occ.size()-1;
    int top=0, bottom=n-1;
    for (int i=(int)pat.size()-1; i>=0; i--) {
        int c = cidx(pat[i]);
        if (c <= 0) return {1,0};
        top    = fm.C[c] + fm.Occ[top][c];
        bottom = fm.C[c] + fm.Occ[bottom+1][c]-1;
        if (top > bottom) return {1,0};
    }
    return {top, bottom};
}

int bwtSearch(const string& ref, const vector<int>& sa,
              const FMIndex& fm, const string& read) {
    const int L = (int)read.size();
    const int N = (int)ref.size();
    int best_pos=-1, best_mm=MAX_MISMATCH+1;
    for (int s = 0; s <= MAX_MISMATCH && best_pos==-1; s++) {
        int seed_offset = s * SEED_LEN;
        if (seed_offset + SEED_LEN > L) break;
        string seed = read.substr(seed_offset, SEED_LEN);
        auto [top, bottom] = backwardSearch(seed, fm);
        if (top > bottom) continue;
        int limit = min(bottom, top+MAX_CAND-1);
        for (int k = top; k <= limit; k++) {
            int cand = sa[k] - seed_offset;
            if (cand < 0 || cand+L > N) continue;
            int mm = 0;
            for (int j = 0; j < L; j++)
                if (ref[cand+j] != read[j])
                    if (++mm >= best_mm) break;
            if (mm < best_mm) { best_mm=mm; best_pos=cand; if(!best_mm) return best_pos; }
        }
    }
    return best_pos;
}

string loadFile(const string& fn) {
    ifstream f(fn);
    if (!f.is_open()) return "";
    string s; getline(f, s); return s;
}

vector<Read> loadReads(const string& fn) {
    ifstream f(fn); vector<Read> reads;
    if (!f.is_open()) return reads;
    reads.reserve(100000);
    string line; getline(f, line);
    while (getline(f, line)) {
        size_t t1=line.find('\t'), t2=line.find('\t',t1+1);
        Read r;
        r.start_pos = stoi(line.substr(t1+1, t2-t1-1));
        r.seq = line.substr(t2+1);
        reads.push_back(r);
    }
    return reads;
}

// 실험 1회 실행 → CSV 한 줄 기록
void runExperiment(
    const string& dataset, const string& snp,
    const string& ref, const vector<int>& sa, const FMIndex& fm,
    const vector<Read>& reads, const string& original,
    ofstream& csv)
{
    int N = (int)ref.size();
    string reconstructed(N, '-');
    int total=0, mapped=0, correct_pos=0;

    clock_t t_start = clock();
    for (auto& r : reads) {
        if ((int)r.seq.size() != 30) continue;
        total++;
        int pos = bwtSearch(ref, sa, fm, r.seq);
        if (pos != -1) {
            mapped++;
            if (pos == r.start_pos) correct_pos++;
            for (int j = 0; j < 30 && pos+j < N; j++)
                reconstructed[pos+j] = r.seq[j];
        }
    }
    clock_t t_end = clock();

    int mismatched = 0;
    for (int i = 0; i < N; i++)
        if (reconstructed[i] != original[i]) mismatched++;

    double total_sec    = (double)(t_end - t_start) / CLOCKS_PER_SEC;
    double memory       = check_memory();
    double mapped_pct   = 100.0 * mapped / total;
    double reconstruct  = 100.0 * (N - mismatched) / N;

    cout << fixed << setprecision(2);
    cout << "[" << dataset << " / SNP " << snp << "]\n";
    cout << "  걸린 시간      : " << total_sec   << " 초\n";
    cout << "  메모리         : " << memory      << " MB\n";
    cout << "  매핑률         : " << mapped_pct  << "%\n";
    cout << "  재구축 일치율  : " << reconstruct << "%\n\n";

    csv << fixed << setprecision(2)
        << "BWT+FM-index," << dataset << "," << snp << ","
        << total_sec << "," << memory << ","
        << mapped_pct << "," << reconstruct << "\n";
}

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif

    // ── CSV 초기화 ─────────────────────────────────────────────
    ofstream csv("results_bwt.csv");
    csv << "algorithm,dataset,snp_rate,total_sec,memory_mb,mapped_pct,reconstruct_pct\n";

    // ── SNP 비율 목록 ──────────────────────────────────────────
    vector<pair<string,string>> snp_list = {
        {"0.1%", "reference_0.1pct.txt"},
        {"1.0%", "reference_1.0pct.txt"},
        {"2.0%", "reference_2.0pct.txt"},
        {"5.0%", "reference_5.0pct.txt"}
    };

    // ── Read 조건 목록 (synthetic) ────────────────────────────
    vector<pair<string,string>> read_list = {
        {"Baseline",  "reads_baseline.txt"},
        {"InDel",     "reads_indel.txt"},
        {"End-Heavy", "reads_end_heavy.txt"}
    };

    string original_synthetic = loadFile("original_synthetic_1M.txt");
    string original_yeast     = loadFile("original_yeast_1M.txt");

    // ── Synthetic 실험 (12가지) ───────────────────────────────
    cout << "===== Synthetic 실험 (12가지) =====\n\n";
    for (auto& [snp_label, ref_file] : snp_list) {
        string ref = loadFile(ref_file);
        if (ref.empty()) {
            cout << "[건너뜀] " << ref_file << " 없음\n\n"; continue;
        }
        cout << "▶ 인덱스 구축 중 (SNP " << snp_label << ")...\n";
        clock_t ti = clock();
        string refS = ref + '$';
        vector<int> sa = buildSA(refS);
        string bwt = buildBWT(refS, sa);
        FMIndex fm(bwt);
        cout << "  인덱스 완료 ("
             << (double)(clock()-ti)/CLOCKS_PER_SEC << "초)\n\n";

        for (auto& [ds_label, reads_file] : read_list) {
            vector<Read> reads = loadReads(reads_file);
            if (reads.empty()) {
                cout << "[건너뜀] " << reads_file << " 없음\n\n"; continue;
            }
            runExperiment(ds_label, snp_label, ref, sa, fm,
                          reads, original_synthetic, csv);
        }
    }

    // ── Yeast 실험 (4가지) ────────────────────────────────────
    if (original_yeast.empty()) {
        cout << "[건너뜀] original_yeast_1M.txt 없음 - yeast 실험 스킵\n";
    } else {
        cout << "===== 빵효모 실험 (4가지) =====\n\n";
        vector<pair<string,string>> yeast_snp_list = {
            {"0.1%", "reference_yeast_0.1pct.txt"},
            {"1.0%", "reference_yeast_1.0pct.txt"},
            {"2.0%", "reference_yeast_2.0pct.txt"},
            {"5.0%", "reference_yeast_5.0pct.txt"}
        };
        vector<Read> yeast_reads = loadReads("reads_yeast.txt");

        for (auto& [snp_label, ref_file] : yeast_snp_list) {
            string ref = loadFile(ref_file);
            if (ref.empty()) {
                cout << "[건너뜀] " << ref_file << " 없음\n\n"; continue;
            }
            cout << "▶ 인덱스 구축 중 (Yeast SNP " << snp_label << ")...\n";
            clock_t ti = clock();
            string refS = ref + '$';
            vector<int> sa = buildSA(refS);
            string bwt = buildBWT(refS, sa);
            FMIndex fm(bwt);
            cout << "  인덱스 완료 ("
                 << (double)(clock()-ti)/CLOCKS_PER_SEC << "초)\n\n";

            runExperiment("빵효모", snp_label, ref, sa, fm,
                          yeast_reads, original_yeast, csv);
        }
    }

    cout << "===================================\n";
    cout << "전체 완료 → results_bwt.csv 저장됨\n";
    cout << "===================================\n";
    return 0;
}
