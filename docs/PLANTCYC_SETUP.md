# PlantCyc/BioCyc API Setup Guide

## 계정 생성 및 인증 설정

PlantCyc API 접근을 위해서는 BioCyc 계정이 필요합니다. 무료 계정 생성 시 1개월 trial이 제공됩니다.

### 1. BioCyc 무료 계정 생성

**계정 생성 페이지**: https://biocyc.org/new-account.shtml

1. 위 링크를 방문하여 계정을 생성합니다
2. 이메일 주소와 비밀번호를 입력합니다
3. 계정 생성 후 이메일 인증을 완료합니다

**혜택**:
- 1개월 무료 trial (모든 BioCyc 데이터베이스 접근 가능)
- 20,000+ 데이터베이스 무제한 접근
- MetaCyc, PlantCyc, SoyCyc 등 식물 대사경로 데이터베이스 포함

### 2. Credential 환경변수 설정

계정 생성 후, 아래와 같이 환경변수를 설정하세요:

```bash
# Linux/Mac
export BIOCYC_EMAIL="your_email@example.com"
export BIOCYC_PASSWORD="your_password"

# 또는 ~/.bashrc 또는 ~/.bash_profile에 추가
echo 'export BIOCYC_EMAIL="your_email@example.com"' >> ~/.bashrc
echo 'export BIOCYC_PASSWORD="your_password"' >> ~/.bashrc
source ~/.bashrc
```

```bash
# Windows PowerShell
$env:BIOCYC_EMAIL = "your_email@example.com"
$env:BIOCYC_PASSWORD = "your_password"

# 또는 시스템 환경변수로 설정 (제어판 > 시스템 > 고급 시스템 설정)
```

### 3. PlantCyc 매핑 실행

환경변수 설정 후 PlantCyc 매핑을 실행합니다:

```bash
cd /data/ethylene
python src/plantcyc_api.py map
```

**예상 소요 시간**: 약 3-5분 (API rate limit: 1 request/sec)

**출력 파일**:
- `data/processed/plantcyc_metabolite_pathways.csv`: 대사체-경로 매핑 결과

### 4. 매핑 결과 확인

```bash
# 매핑 통계 확인
head -20 data/processed/plantcyc_metabolite_pathways.csv

# 매핑 성공률 확인
python -c "
import pandas as pd
df = pd.read_csv('data/processed/plantcyc_metabolite_pathways.csv')
total = len(df)
mapped = df['PlantCyc_Pathway_ID'].notna().sum()
print(f'Total: {total}, Mapped: {mapped}, Rate: {100*mapped/total:.1f}%')
"
```

### 5. PlantCyc Pathway Enrichment 분석

매핑 완료 후 pathway enrichment 분석을 수행합니다:

```bash
python src/plantcyc_pathway_enrichment.py
```

**출력 파일**:
- `results/plantcyc_pathway_enrichment.csv`: Enriched pathway 목록

---

## Troubleshooting

### 문제: "Failed to parse JSON response"

**원인**: 인증이 실패하거나 세션이 만료됨

**해결**:
1. 환경변수가 올바르게 설정되었는지 확인
   ```bash
   echo $BIOCYC_EMAIL
   echo $BIOCYC_PASSWORD
   ```
2. 계정 정보가 정확한지 확인 (웹사이트에서 로그인 테스트)
3. 스크립트를 다시 실행

### 문제: "Request failed: 404"

**원인**: ORGID가 잘못되었거나 데이터베이스가 존재하지 않음

**해결**:
- **META (MetaCyc)** 사용 (가장 포괄적인 데이터베이스)
- ARA (AraCyc) 또는 다른 organism-specific database 시도

### 문제: 매핑률이 매우 낮음

**원인**:
- 대사체 이름이 BioCyc에 등록되지 않음
- 식물 특이적 대사체 (예: malonyl conjugates)가 MetaCyc에 없을 수 있음

**해결**:
1. ChEBI ID 기반 매핑 시도 (향후 구현 예정)
2. KEGG ID 기반 cross-reference
3. 수동으로 주요 대사체 (Daidzein, Formononetin 등) 매핑

---

## Alternative: KEGG Pathway Analysis (권장)

PlantCyc 접근이 어려운 경우, 이미 완료된 KEGG 분석 결과를 사용할 수 있습니다:

**KEGG 결과 파일**: `results/table1_metabolomics_real.csv`

**주요 발견**:
- **map01110** (Biosynthesis of secondary metabolites): **P = 0.030** (유의미!)
- 5개의 significant 대사체가 모두 이 경로에 속함
- Isoflavonoid 경로 활성화 확인

KEGG 분석만으로도 논문 작성이 가능하며, PlantCyc는 추가 validation으로 사용할 수 있습니다.

---

## References

- BioCyc Web Services: https://biocyc.org/web-services.shtml
- SoyCyc Database: https://plantcyc.org/databases/soycyc/
- Plant Metabolic Network: https://plantcyc.org/
- MetaCyc: https://metacyc.org/
