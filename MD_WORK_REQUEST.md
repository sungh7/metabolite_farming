# MD 시뮬레이션 업무 지시서

**프로젝트**: 에틸렌 유도 이소플라보노이드 생합성 연구
**작성일**: 2026년 1월 15일
**요청자**:
**수행자**: MD 담당자

---

## 1. 배경 및 목적

에틸렌 처리된 콩(Glycine max)의 멀티오믹스 분석 결과, 말로닐화/아세틸화 이소플라본 접합체가 4,000배 이상 증가함을 확인하였습니다. 이들 대사체와 생합성 효소 간의 결합을 분자 도킹으로 예측한 결과, **문헌에 보고되지 않은 신규 결합 조합**을 발견하였습니다.

**목적**: MD 시뮬레이션을 통해 도킹 결과의 결합 안정성 및 결합 자유 에너지 검증

---

## 2. 우선순위별 MD 시뮬레이션 대상

### 최우선 (Priority 1) - 완전 신규

| 순위 | 리간드 | 수용체 | PDB | 도킹 점수 | 신규성 |
|------|--------|--------|-----|-----------|--------|
| **1** | **6''-O-Malonyldaidzin** | β-Glucosidase | 6YN7 | **-9.165 kcal/mol** | ★★★ |
| **2** | **6''-O-Acetyldaidzin** | 2-HIS (IFS) | 8E83 | **-8.863 kcal/mol** | ★★★ |
| **3** | **6''-O-Acetylgenistin** | 2-HIS (IFS) | 8E83 | **-7.768 kcal/mol** | ★★★ |
| **4** | **6''-O-Malonylgenistin** | 2-HIS (IFS) | 8E83 | **-7.485 kcal/mol** | ★★★ |

### 우선 (Priority 2) - 부분 신규

| 순위 | 리간드 | 수용체 | PDB | 도킹 점수 | 비고 |
|------|--------|--------|-----|-----------|------|
| 5 | Daidzin | β-Glucosidase | 6YN7 | -8.496 kcal/mol | Km 데이터 있음 |
| 6 | Genistin | β-Glucosidase | 6YN7 | -8.500 kcal/mol | Km 데이터 있음 |
| 7 | Daidzin | 2-HID | 8EA1 | -7.357 kcal/mol | 기질 |

---

## 3. 구조 파일 위치

### 수용체 (PDB)
```
/data/ethylene/data/structures/pdb/
├── 6YN7.pdb    # β-Glucosidase (1.27 MB) - Alicyclobacillus herbarius
├── 8E83.pdb    # 2-HIS/IFS homolog (1.32 MB) - Medicago truncatula, 2.0 Å
├── 8EA1.pdb    # 2-HID (445 KB) - Pueraria lobata, 2.4 Å
└── 1EYQ.pdb    # CHI (345 KB) - Medicago sativa, 1.85 Å
```

### 리간드 (SDF)
```
/data/ethylene/data/structures/ligands/
├── 6-O-Malonyldaidzin_CID5318574.sdf
├── 6-O-Malonylgenistin_CID5318568.sdf
├── 6-O-Acetyldaidzin_CID14034712.sdf
├── 6-O-Acetylgenistin_CID5320413.sdf
├── Daidzin_CID107971.sdf
├── Genistin_CID5281377.sdf
├── Daidzein_CID5281708.sdf
└── Genistein_CID5280961.sdf
```

### 도킹 결과 (초기 포즈)
```
/data/ethylene/results/docking/all_candidates/
├── 6YN7_6-O-Malonyldaidzin_CI/docked.pdbqt   # 최우선
├── 8E83_6-O-Acetyldaidzin_CID1403/docked.pdbqt
├── 8E83_6-O-Acetylgenistin_CID532/docked.pdbqt
└── [기타 조합별 폴더]/
```

---

## 4. 요청 작업

### 4.1 필수 작업

1. **100 ns Production MD** (최소 3개 시스템)
   - Priority 1 중 상위 3개 조합
   - AMBER 또는 CHARMM force field 사용
   - 리간드 파라미터: GAFF2 또는 CGenFF

2. **결합 자유 에너지 계산**
   - MM-PBSA 또는 MM-GBSA
   - 마지막 50 ns 구간 사용
   - 잔기별 분해 분석 포함

3. **결합 안정성 분석**
   - Ligand RMSD (수렴 여부)
   - Protein-Ligand 수소결합 개수
   - 결합 부위 RMSF

### 4.2 선택 작업

- 200 ns 연장 시뮬레이션 (안정성 미확인 시)
- Alanine scanning (핵심 잔기 검증)
- 비교 시뮬레이션 (Daidzin vs Malonyldaidzin)

---

## 5. 예상 산출물

| 산출물 | 형식 | 설명 |
|--------|------|------|
| RMSD plot | PNG/PDF | 리간드 RMSD 시간 변화 |
| Binding free energy | 표 | ΔG_bind (kcal/mol) ± 오차 |
| Key residues | 목록 | 상호작용 기여도 상위 10개 |
| H-bond analysis | 그래프 | 수소결합 개수 시간 변화 |
| Representative snapshot | PDB | 결합 상태 대표 구조 |

---

## 6. 참고 정보

### 6.1 기존 문헌 데이터 (비교용)

| 시스템 | 결합 에너지 | 출처 |
|--------|-------------|------|
| Naringenin + CHI | ~ -7.5 kcal/mol | PDB 1EYQ 결정 구조 |
| Daidzin + GmIMaT | Km = 36.28 μM | Ahmad et al. 2017 |
| Genistin + GmIMaT | Km = 23.04 μM | Ahmad et al. 2017 |

### 6.2 핵심 결합 잔기 (도킹 분석 결과)

**6''-O-Acetyldaidzin + 8E83 (IFS)**:
- 수소결합: ASP230, TYR212, LYS223, THR215
- 소수성: ILE45, LEU48, PHE224, PHE237
- 예상 결합 에너지: -8.9 kcal/mol

### 6.3 MDP 파일 위치
```
/data/ethylene/md_protocol/mdp/
├── em.mdp          # Energy minimization
├── nvt.mdp         # NVT equilibration (100 ps)
├── npt.mdp         # NPT equilibration (100 ps)
└── md_100ns.mdp    # Production MD (100 ns)
```

---

## 7. 일정

| 단계 | 예상 소요 | 비고 |
|------|-----------|------|
| 시스템 준비 | 1일 | 리간드 파라미터화 포함 |
| Equilibration | 0.5일 | EM + NVT + NPT |
| Production MD | 3-5일 | GPU 사용 시 |
| 분석 | 1-2일 | MM-PBSA 포함 |
| **총 예상** | **5-8일** | 시스템 3개 기준 |

---

## 8. 연락처

- **프로젝트 관련 문의**:
- **데이터 위치**: `/data/ethylene/`
- **GitHub**: https://github.com/sungh7/metabolite_farming

---

## 9. 체크리스트

- [ ] 수용체 PDB 파일 확인
- [ ] 리간드 SDF 파일 확인
- [ ] 도킹 포즈 (PDBQT) 확인
- [ ] Force field 선택
- [ ] 리간드 파라미터 생성
- [ ] MD 시스템 구축
- [ ] Production MD 실행
- [ ] 분석 완료
- [ ] 결과 보고

---

**작성**: Claude Code (AI Assistant)
**검토**:
