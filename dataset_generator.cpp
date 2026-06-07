// =============================================================
//  Short-read re-sequencing 프로젝트 - 역할3 (이동건)
//  데이터셋 일괄 생성기 (SNP 비율별 reference 자동 생성)
//
//  [목적]
//   main(reference_builder) 을 SNP 비율마다 다시 빌드/실행할 필요 없이,
//   한 번 실행으로 SNP 4종(0.1/1/2/5%)의 reference 와 공용 reads 를 모두 만든다.
//   -> hash_greedy / hash_dp 가 한 번 실행으로 전체 매트릭스를 돌릴 수 있게 한다.
//
//  [중요: 기존 generator 와 동일 산출물]
//   reference_builder.cpp / trap_generator.cpp / read_generator.cpp 의 로직을
//   그대로 복제(시드/RNG 호출 순서 동일)했다. 따라서 SEED=42 기준으로
//   기존 산출물과 byte 단위로 동일하다(= brute-force 표와 공정 비교 가능).
//
//  [생성 파일]
//   원본(1회):    original_synthetic_1M.txt, original_yeast_1M.txt(sequence.fasta 있을 때)
//   공용 reads:   reads_baseline.txt / reads_indel.txt / reads_end_heavy.txt / reads_yeast.txt
//                 (reads 는 무변이 원본에서 잘라내므로 SNP 와 무관 -> 1벌만 생성)
//   SNP별 ref:    reference_synthetic_snp{1,10,20,50}.txt (+ snp_list_synthetic_snp{L}.txt)
//                 reference_yeast_snp{1,10,20,50}.txt     (+ snp_list_yeast_snp{L}.txt)
//                 라벨 규칙(per-mille): snp1=0.1%, snp10=1%, snp20=2%, snp50=5%
//
//  [실행 위치] CLion 은 cmake-build-debug/ 에서 실행하며 sequence.fasta 도 거기에 있다.
// =============================================================

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <array>

using namespace std;

// ──────────────── 파라미터 (기존 generator 와 동일) ────────────────
static const size_t   N           = 1000000;  // 원본 길이
static const unsigned SEED        = 42;        // 고정 시드
static const int      READ_LENGTH = 30;        // read 길이 L
static const int      NUM_READS   = 100000;    // read 수 M
static const int      SEED_BPE    = 20;        // End-Heavy seed 구간 (앞 20bp)
static const string   INPUT_FASTA = "sequence.fasta";

static const array<char, 4> BASES = {'A', 'C', 'G', 'T'};

// SNP 라벨/비율 테이블
struct SnpLevel { string label; double rate; };
static const vector<SnpLevel> SNP_LEVELS = {
    {"1",  0.001},  // 0.1%
    {"10", 0.01},   // 1.0%
    {"20", 0.02},   // 2.0%
    {"50", 0.05},   // 5.0%
};

// ──────────────── reference_builder.cpp 복제 ────────────────
static string generateSynthetic(size_t need, mt19937& rng) {
    uniform_int_distribution<int> pick(0, 3);
    string seq;
    seq.reserve(need);
    for (size_t i = 0; i < need; ++i) seq.push_back(BASES[pick(rng)]);
    return seq;
}

static char mutateBase(char base, mt19937& rng) {
    int idx;
    switch (base) {
        case 'A': idx = 0; break;
        case 'C': idx = 1; break;
        case 'G': idx = 2; break;
        case 'T': idx = 3; break;
        default:  return base;
    }
    uniform_int_distribution<int> offset(1, 3);
    return BASES[(idx + offset(rng)) & 3];
}

static string extractOriginal(const string& fastaPath, size_t need) {
    ifstream in(fastaPath);
    if (!in.is_open()) return "";
    string seq;
    seq.reserve(need);
    string line;
    while (getline(in, line) && seq.size() < need) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) {
            if (c == 'a' || c == 'A') c = 'A';
            else if (c == 'c' || c == 'C') c = 'C';
            else if (c == 'g' || c == 'G') c = 'G';
            else if (c == 't' || c == 'T') c = 'T';
            else continue;
            seq.push_back(c);
            if (seq.size() >= need) break;
        }
    }
    return seq;
}

// original 에 SNP 를 심어 reference / snp_list 저장 (mutateAndSave 복제, 파일명 파라미터화)
static bool mutateAndSave(const string& original, double snp_rate, unsigned seed,
                          const string& out_reference, const string& out_snp_list,
                          const string& label) {
    string reference = original;
    mt19937 rng(seed);
    uniform_real_distribution<double> prob(0.0, 1.0);

    ofstream snpOut(out_snp_list);
    if (!snpOut.is_open()) { cerr << "[오류] " << out_snp_list << " 생성 실패.\n"; return false; }
    snpOut << "position\toriginal\treference\n";

    size_t snpCount = 0;
    for (size_t i = 0; i < reference.size(); ++i) {
        if (prob(rng) < snp_rate) {
            char before = reference[i];
            char after  = mutateBase(before, rng);
            reference[i] = after;
            snpOut << i << '\t' << before << '\t' << after << '\n';
            ++snpCount;
        }
    }

    ofstream out(out_reference);
    if (!out.is_open()) { cerr << "[오류] " << out_reference << " 생성 실패.\n"; return false; }
    out << reference;

    double actualRate = (double)snpCount / reference.size() * 100.0;
    cout << "       [" << label << "] " << out_reference
         << " (SNP " << snpCount << "개, 실제 " << actualRate << "%)\n";
    return true;
}

// ──────────────── trap_generator.cpp 복제 (synthetic reads) ────────────────
static char randomBase(mt19937& rng) {
    return BASES[uniform_int_distribution<int>(0, 3)(rng)];
}
static char mutate(char c, mt19937& rng) {
    int idx = (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : 3;
    return BASES[(idx + uniform_int_distribution<int>(1, 3)(rng)) & 3];
}

static void genBaseline(const string& orig, mt19937& rng, const string& out) {
    ofstream fout(out);
    fout << "read_id\tstart_pos\tsequence\n";
    uniform_int_distribution<int> pos(0, (int)orig.size() - READ_LENGTH);
    for (int i = 0; i < NUM_READS; i++) {
        int start = pos(rng);
        fout << i << '\t' << start << '\t' << orig.substr(start, READ_LENGTH) << '\n';
    }
    cout << "       " << out << " (에러 없음)\n";
}

static void genInDel(const string& orig, mt19937& rng, const string& out) {
    ofstream fout(out);
    fout << "read_id\tstart_pos\tsequence\n";
    uniform_int_distribution<int> pos(0, (int)orig.size() - READ_LENGTH);
    uniform_int_distribution<int> indelCount(1, 2);
    uniform_int_distribution<int> indelType(0, 1);
    for (int i = 0; i < NUM_READS; i++) {
        int start = pos(rng);
        string seg = orig.substr(start, READ_LENGTH);
        int n = indelCount(rng);
        uniform_int_distribution<int> indelPos(0, (int)seg.size() - 1);
        for (int k = 0; k < n; k++) {
            int p = indelPos(rng) % (int)seg.size();
            if (indelType(rng) == 0) seg.insert(seg.begin() + p, randomBase(rng));
            else if (seg.size() > 1)  seg.erase(seg.begin() + p);
        }
        fout << i << '\t' << start << '\t' << seg << '\n';
    }
    cout << "       " << out << " (InDel 1~2개, 치환 없음)\n";
}

static void genEndHeavy(const string& orig, mt19937& rng, const string& out) {
    ofstream fout(out);
    fout << "read_id\tstart_pos\tsequence\n";
    uniform_int_distribution<int> pos(0, (int)orig.size() - READ_LENGTH);
    uniform_int_distribution<int> tailPos(SEED_BPE, READ_LENGTH - 1);
    uniform_int_distribution<int> errType(0, 2);
    for (int i = 0; i < NUM_READS; i++) {
        int start = pos(rng);
        string seg = orig.substr(start, READ_LENGTH);
        vector<int> usedPos;
        int injected = 0, attempts = 0;
        while (injected < 2 && attempts < 20) {
            attempts++;
            int p = tailPos(rng);
            if (p >= (int)seg.size()) continue;
            bool dup = false;
            for (int u : usedPos) if (u == p) { dup = true; break; }
            if (dup) continue;
            int et = errType(rng);
            if (et == 0)       seg[p] = mutate(seg[p], rng);
            else if (et == 1)  seg.insert(seg.begin() + p, randomBase(rng));
            else {
                if (seg.size() > SEED_BPE + 1) seg.erase(seg.begin() + p);
                else seg[p] = mutate(seg[p], rng);
            }
            usedPos.push_back(p);
            injected++;
        }
        fout << i << '\t' << start << '\t' << seg << '\n';
    }
    cout << "       " << out << " (앞 20bp 완벽, 뒤 10bp 에러 2개)\n";
}

// ──────────────── read_generator.cpp 복제 (clean reads, yeast 용) ────────────────
static void genCleanReads(const string& orig, const string& out, const string& label) {
    ofstream fout(out);
    fout << "read_id\tstart_pos\tsequence\n";
    mt19937 rng(SEED);
    uniform_int_distribution<int> startDist(0, (int)orig.size() - READ_LENGTH);
    for (int i = 0; i < NUM_READS; ++i) {
        int start = startDist(rng);
        fout << i << '\t' << start << '\t' << orig.substr(start, READ_LENGTH) << '\n';
    }
    cout << "       [" << label << "] " << out << " (clean reads)\n";
}

static void writeOriginal(const string& seq, const string& out) {
    ofstream o(out);
    o << seq;
}

// 원본 1개에 대해: 공용 reads + SNP별 reference 생성
static void buildAll(const string& original, bool is_synthetic,
                     const string& ref_prefix, const string& snplist_prefix) {
    // SNP별 reference
    for (const auto& s : SNP_LEVELS) {
        string ref = ref_prefix + "_snp" + s.label + ".txt";
        string snplist = snplist_prefix + "_snp" + s.label + ".txt";
        mutateAndSave(original, s.rate, SEED, ref, snplist, ref_prefix);
    }
}

int main() {
    cout << "=== 데이터셋 일괄 생성 시작 (SEED=" << SEED << ") ===\n";

    // ── A. 인공 서열 ──────────────────────────────────────────
    cout << "[A] 인공 원본 생성 (mt19937, ATCG 균등, " << N << "bp)\n";
    mt19937 synth_rng(SEED);
    string synthetic = generateSynthetic(N, synth_rng);
    writeOriginal(synthetic, "original_synthetic_1M.txt");
    cout << "       original_synthetic_1M.txt 생성 완료\n";

    cout << "    [공용 reads] synthetic 3종 생성\n";
    {
        mt19937 rng1(SEED), rng2(SEED + 1), rng3(SEED + 2);  // trap_generator 와 동일
        genBaseline (synthetic, rng1, "reads_baseline.txt");
        genInDel    (synthetic, rng2, "reads_indel.txt");
        genEndHeavy (synthetic, rng3, "reads_end_heavy.txt");
    }
    cout << "    [SNP별 reference] reference_synthetic_snp{1,10,20,50}.txt\n";
    buildAll(synthetic, true, "reference_synthetic", "snp_list_synthetic");

    // ── B. 빵효모 (sequence.fasta 있을 때만) ──────────────────
    ifstream test(INPUT_FASTA);
    if (test.is_open()) {
        test.close();
        cout << "\n[B] " << INPUT_FASTA << " 발견 -> 빵효모 처리\n";
        string yeast = extractOriginal(INPUT_FASTA, N);
        if (yeast.size() < N) {
            cerr << "[경고] 빵효모 추출 부족 (" << yeast.size() << "/" << N << ") -> 건너뜀.\n";
        } else {
            writeOriginal(yeast, "original_yeast_1M.txt");
            cout << "       original_yeast_1M.txt 생성 완료 (" << yeast.size() << "bp)\n";
            cout << "    [공용 reads] yeast clean reads 생성\n";
            genCleanReads(yeast, "reads_yeast.txt", "yeast");
            cout << "    [SNP별 reference] reference_yeast_snp{1,10,20,50}.txt\n";
            buildAll(yeast, false, "reference_yeast", "snp_list_yeast");
        }
    } else {
        cout << "\n[알림] " << INPUT_FASTA << " 없음 -> 빵효모 건너뜀 (synthetic 만).\n";
    }

    cout << "\n=== 생성 완료 ===\n";
    return 0;
}
