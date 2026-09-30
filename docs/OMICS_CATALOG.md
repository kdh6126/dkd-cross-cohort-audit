# 공개 자료 카탈로그

무엇을 어디서 받았고, 참여자 한 명이 어떤 오믹스 층을 가지고 있으며,
그래서 **같은 참여자에서 두 층을 잰 비교가 성립하는가**를 조회한다.

## 만들기

    python scripts/build_omics_catalog.py      # db/omics_catalog.sqlite

## 브라우저에서 보기

    python -m streamlit run scripts/omics_catalog_app.py

http://localhost:8501 이 열린다.

## 표와 뷰

| 이름 | 내용 |
|---|---|
| `dataset` | 접근번호, 저장소, 오믹스, 조직, 사례/대조 수, 조달 경로 |
| `subject_layer` | 참여자별 보유 층(전사체·단백체·대사체), 질환 분류, 사례/대조 |
| `v_multiomics_feasibility` | 층 쌍마다 사례와 대조가 각각 몇 명인가 |
| `v_layer_by_class` | 질환 분류별 층 보유 인원과 다층 보유자 수 |
| `v_dataset_by_omics` | 오믹스별 자료 수와 인원 |

## 판정 규칙

다중 오믹스가 성립하려면 **한 층 쌍에 대해 사례와 대조 양쪽에 사람이 있어야** 한다.
한쪽이 0명이면 짝은 존재해도 비교가 안 된다. 이 프로젝트가 두 번 틀린 지점이라
사람이 세지 않고 뷰가 세게 두었다.

## 주의

참여자 식별자는 KPMP 공개 아틀라스의 비식별 번호다. 그래도 공개본(`release/`)에는
넣지 않는다. 논문이 쓰는 것은 질환별 집계뿐이므로 개별 행을 배포할 이유가 없다.
