# Short-read Re-sequencing 프로젝트

Reference 기반 (seed-and-extend, 비둘기집 원리) short-read 매칭을 위한 **공통 데이터 생성 파이프라인**.

팀원 4명이 각자 다른 알고리즘(Hash+Greedy / RB-tree+DP / Minimizer+SA / BWT+SA)으로
**같은 입력 데이터**에서 시간·정확도를 비교할 수 있도록, 깨끗한 원본/레퍼런스/리드를 생성한다.

---

## 빠른 실행 가이드

1. CLion 에서 빌드 후 실행 (또는 터미널에서 `./2026_algorithm_project`)
2. 끝. `cmake-build-debug/` 폴더에 결과 파일들이 자동 생성된다.
3. **빵효모 데이터를 비교 실험에 추가하고 싶다면** `cmake-build-debug/` 폴더에 `sequence.fasta` 파일을 미리 넣어두면 자동으로 함께 처리된다.

> 인공 서열은 mt19937 난수로 매번 동일하게 생성되므로 별도 입력 파일은 불필요.

---

## 파일 구조

### `main.cpp` — 파이프라인 진입점

- 역할: `build_reference()` → `generate_reads()` 순서로 호출
- **SNP 비율은 여기 한 줄만 바꾸면 된다**
  ```cpp
  static const double SNP_RATE = 0.001;  // 0.001=0.1% / 0.005=0.5% / 0.01=1%
  ```

### `reference_builder.cpp` — 1단계 (역할1, 이동건)

`build_reference(double snp_rate)` 호출 시 두 종류의 원본을 만들고 각각에 SNP 를 심는다.

**실행 결과로 생성되는 파일**

| 출력 파일                    | 종류                  | 내용                                            | 생성 조건                |
| ---------------------------- | --------------------- | ----------------------------------------------- | ------------------------ |
| `original_synthetic_1M.txt`  | 인공 원본             | mt19937 으로 ATCG 25% 균등 랜덤 1,000,000 염기  | **항상**                 |
| `reference_synthetic.txt`    | 인공 reference (+SNP) | 위 원본에 `snp_rate` 확률로 SNP 삽입            | **항상**                 |
| `snp_list_synthetic.txt`     | 인공 정답표           | `position \t original \t reference` 형식        | **항상**                 |
| `original_yeast_1M.txt`      | 빵효모 원본           | `sequence.fasta` 에서 ATCG 만 추출한 100만 염기 | `sequence.fasta` 있을 때 |
| `reference_yeast.txt`        | 빵효모 reference      | 위 원본에 동일 `snp_rate` 로 SNP 삽입           | `sequence.fasta` 있을 때 |
| `snp_list_yeast.txt`         | 빵효모 정답표         | 동일 형식                                       | `sequence.fasta` 있을 때 |

- 두 원본에 **같은 SEED(=42)** 로 SNP 를 심으므로 변이 위치가 동일 → 인공 vs 빵효모 직접 비교가 공정
- **인공 서열을 baseline 으로 쓰는 이유**: 실제 유전체(빵효모)는 repeat 구역과 GC content 편향이 심해 알고리즘의 순수 성능 비교에 부적합
- 난수: `std::mt19937(SEED=42)` 고정 → 팀원 모두 동일한 결과

### `read_generator.cpp` — 2단계 (역할2, 김세훈)

`generate_reads()` 호출 시 위에서 만들어진 원본 파일을 자동 감지해 두 세트의 read 를 생성한다.

**실행 결과로 생성되는 파일**

| 출력 파일             | 입력 원본                   | 내용                  | 생성 조건             |
| --------------------- | --------------------------- | --------------------- | --------------------- |
| `reads_synthetic.txt` | `original_synthetic_1M.txt` | 30bp read 100,000개   | **항상**              |
| `reads_yeast.txt`     | `original_yeast_1M.txt`     | 30bp read 100,000개   | 빵효모 원본 있을 때만 |

- read 자체에는 에러를 넣지 않음 (변이는 reference 쪽에만 존재)
- 형식: `read_id \t start_pos \t sequence` (start_pos 는 해당 original 기준 0-based 정답 위치)
- 난수: `std::mt19937(SEED=42)`, `uniform_int_distribution` 사용
  → `rand()` 의 짧은 주기/패턴 문제 없음

### `trap_generator.cpp` — Trap Read 생성기 (역할2, 김세훈)

12개 실험 조합을 위한 **3가지 조건의 read 데이터셋** 생성. 별도 실행파일로 독립 실행.

```bash
g++ -O2 -std=c++17 -o trap_generator trap_generator.cpp
./trap_generator
```

> `original_synthetic_1M.txt` 가 같은 폴더에 있어야 함 (main pipeline 먼저 실행)

**실행 결과로 생성되는 파일**

| 출력 파일              | 조건                              | 실험 목적                                      |
| ---------------------- | --------------------------------- | ---------------------------------------------- |
| `reads_baseline.txt`   | 에러 없음 (깨끗한 read)           | 순수 탐색 속도 비교 기준선                     |
| `reads_indel.txt`      | InDel 1~2개 주입, 치환 없음       | 길이 변화 발생 시 brute-force 격파 검증        |
| `reads_end_heavy.txt`  | 앞 20bp 완벽 + 뒤 10bp 에러 2개  | Seed-and-Extend 전략 타당성 증명               |

- **SNP 비율별 reference 4종 (이동건) × read 조건 3종 = 12개 실험 조합**
- InDel은 read에 주입 (시퀀서 기계 오류 시뮬레이션), SNP는 reference에만 존재
- SEED 독립 적용 (SEED / SEED+1 / SEED+2) → 각 조건 재현 가능

### `CMakeLists.txt`

- `main.cpp` + `reference_builder.cpp` + `read_generator.cpp` → 실행파일 `2026_algorithm_project`
- `trap_generator.cpp` → 별도 실행파일 `trap_generator`

---

## 파라미터 조정 위치

| 항목        | 위치                                  | 기본값                    |
| ----------- | ------------------------------------- | ------------------------- |
| SNP 비율    | `main.cpp` 의 `SNP_RATE`              | 0.001 (0.1%)              |
| 원본 길이 N | `reference_builder.cpp` 의 `N`        | 1,000,000                 |
| Read 길이 L | `read_generator.cpp` 의 `READ_LENGTH` | 30                        |
| Read 개수 M | `read_generator.cpp` 의 `NUM_READS`   | 100,000                   |
| 난수 시드   | 두 cpp 의 `SEED`                      | 42 (반드시 동일하게 유지) |

---

## 결과 파일 포맷

### `original_*_1M.txt`, `reference_*.txt`

- 헤더 없음. ATCG 만으로 구성된 1,000,000 글자 한 줄

### `snp_list_*.txt`

```
position    original    reference
13          A           C
2602        T           A
...
```

### `reads_*.txt`

```
read_id    start_pos    sequence
0          488213       ACGTACGTACGTACGTACGTACGTACGTAC
1          991122       TGCATGCATGCATGCATGCATGCATGCATG
...
```

---

## 다음 단계 (팀원 3·4)

- **기본 실험 입력**: `reference_synthetic.txt` + `reads_synthetic.txt` + 정답 비교용 `snp_list_synthetic.txt`
- **비교 실험 입력** (반복서열·GC 편향 효과 확인): `*_yeast.txt` 세트

`reads_*.txt` 의 각 read 가 해당 `reference_*.txt` 게놈의 어느 위치에서 왔는지 찾는 것이 본 과제의 핵심.
정답이 `reads_*.txt` 의 `start_pos` 컬럼에 있으므로 정확도(recall/precision) 측정 가능.
SNP 가 심긴 위치에 read 가 걸리면 mismatch 1~3개가 발생하므로 비둘기집 원리(L/3 perfect seed) 가 동작한다.
