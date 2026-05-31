# Short-read Re-sequencing 프로젝트

Reference 기반 (seed-and-extend, 비둘기집 원리) short-read 매칭을 위한 **공통 데이터 생성 파이프라인**.

팀원 4명이 각자 다른 알고리즘(Hash+Greedy / RB-tree+DP / Minimizer+SA / BWT+SA)으로
**같은 입력 데이터**에서 시간·정확도를 비교할 수 있도록, 깨끗한 원본/레퍼런스/리드를 생성한다.

---

## 빠른 실행 가이드

1. **빌드 결과 폴더(`cmake-build-debug/`)에 원본 FASTA 파일을 넣는다.**
   - 파일명은 `sequence.fasta` 로 둘 것
   - 예: 빵효모(S. cerevisiae) 1번 염색체 1Mb 분량
2. CLion 에서 빌드 후 실행 (또는 터미널에서 `./2026_algorithm_project`)
3. 끝. 같은 폴더에 4개의 결과 파일이 생성된다.

```
cmake-build-debug/
├── sequence.fasta            <-- 사용자가 직접 넣어야 하는 유일한 파일
├── original_1M.txt           (자동 생성) 원본 100만 염기
├── reference_genome.txt      (자동 생성) SNP 삽입본
├── snp_list.txt              (자동 생성) 변이 위치 정답표
└── reads.txt                 (자동 생성) 길이 30bp read 10만개
```

---

## 파일 구조

### `main.cpp` — 파이프라인 진입점

- 역할: `build_reference()` → `generate_reads()` 순서로 호출만 함
- **SNP 비율은 여기 한 줄만 바꾸면 된다**
  ```cpp
  static const double SNP_RATE = 0.001;  // 0.001=0.1% / 0.005=0.5% / 0.01=1%
  ```

### `reference_builder.cpp` — 1단계 (역할1, 이동건)

- `build_reference(double snp_rate)` 함수 하나로 두 가지를 수행
  1. `sequence.fasta` 에서 ATCG 만 추출해 정확히 1,000,000 염기의 `original_1M.txt` 생성
  2. 원본을 복사한 뒤 `snp_rate` 확률로 각 위치를 다른 염기로 치환 → `reference_genome.txt`
  3. 변이 위치를 `snp_list.txt` 에 `position \t original \t reference` 형식으로 기록
- 난수: `std::mt19937(SEED=42)` 고정 → 팀원 모두 동일 결과

### `read_generator.cpp` — 2단계 (역할2, 김세훈)

- `generate_reads()` 함수
- `original_1M.txt` 에서 랜덤 위치를 잡아 길이 **30bp** read **100,000개** 추출
- read 자체에는 에러를 넣지 않음 (변이는 reference 쪽에만 존재)
- `reads.txt` 형식: `read_id \t start_pos \t sequence`
  - `start_pos` 는 original 기준 0-based 정답 위치 → 매칭 알고리즘 정확도 평가용
- 난수: `std::mt19937(SEED=42)`, `uniform_int_distribution` 사용
  → `rand()` 의 짧은 주기 / 패턴 문제 없음

### `CMakeLists.txt`

- 세 cpp 를 하나의 실행파일 `2026_algorithm_project` 로 묶어 빌드

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

### `original_1M.txt`, `reference_genome.txt`

- 헤더 없음. ATCG 만으로 구성된 1,000,000 글자 한 줄

### `snp_list.txt`

```
position    original    reference
13          A           C
2602        T           A
...
```

### `reads.txt`

```
read_id    start_pos    sequence
0          488213       ACGTACGTACGTACGTACGTACGTACGTAC
1          991122       TGCATGCATGCATGCATGCATGCATGCATG
...
```

---

## 다음 단계 (팀원 3·4)

`reads.txt` 의 각 read 가 `reference_genome.txt` 의 어느 위치에서 왔는지 찾는 것이 본 과제의 핵심.
정답은 `reads.txt` 의 `start_pos` 컬럼에 있으므로 정확도(recall/precision) 측정이 가능하며,
SNP 가 심긴 위치에 read 가 걸리면 mismatch 1~3개가 발생하므로 비둘기집 원리(L/3 perfect seed) 가 동작한다.
