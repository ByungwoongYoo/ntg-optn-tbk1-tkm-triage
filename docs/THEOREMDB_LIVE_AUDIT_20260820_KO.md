# TheoremDB P2624 실시간 공개상태 감사 — 2026-08-20

## 결론

P2624의 공식 상태는 아직 `Open`이며, 공식 상태문은 여전히 `13 <= M <= 14`이다. 다만 현재 문제의 Work 탭에는 초기 packet 기록 5개 외에 최근 기여 12개가 붙어 있고, 그중 제3자 Jordan Boisclair가 별도 정적 인증 패키지를 독립 재생하여 `M=13`을 지지한 R6088–R6090이 있다. 이 기록들은 아직 `Needs packet proposal` 단계이며 공식 상태문에는 반영되지 않았다.

## 확인한 공식 페이지

- 문제: https://www.theoremdb.org/statements/b2-two-set-z100/
- 제3자 시도: https://www.theoremdb.org/record/?kind=research&ref=R6088
- 제3자 주장: https://www.theoremdb.org/record/?kind=research&ref=R6089
- 제3자 아티팩트: https://www.theoremdb.org/record/?kind=research&ref=R6090
- 사용자 계정의 v2 시도: https://www.theoremdb.org/record/?kind=research&ref=R5662
- 사용자 계정의 v2 주장: https://www.theoremdb.org/record/?kind=research&ref=R5663
- 사용자 계정의 v2 아티팩트: https://www.theoremdb.org/record/?kind=research&ref=R5664

## 공식 상태와 최근 기여의 차이

| 층위 | 현재 표시 | 의미 |
|---|---|---|
| 공식 Status 탭 | `Open`; `13 <= M <= 14`; 14-set 존재 여부 미해결 | 검토·packet 반영이 끝난 공식 상태 |
| 초기 packet 기록 | 5개 | R44–R47 및 중복 표시 1개 |
| 최근 기여 | 12개 | packet 이후 붙은 live work이며 공식 상태와 동일하지 않음 |
| Jordan 기록 R6088–R6090 | 완료/지지/사용 가능, 모두 `Needs packet proposal` | 강한 독립 계산 감사지만 아직 packet 검토 전 |

`Established`인 R47이 증명한 것은 정확한 최댓값 13이 아니라 `13 <= M <= 14`이다.

## 사용자 계정의 실제 B2 기록

현재 `@statement` 계정으로 확인되는 B2 관련 기록은 다음 여섯 개다.

| ID | 종류 | 증거 등급 | 저장 상태 | 핵심 |
|---|---|---|---|---|
| R5659 | attempt | computational | Partial | 최초 raw-evidence 보고 |
| R5660 | claim | computational | Supported | 13-set 및 보고된 무-witness 범위의 제한적 주장 |
| R5661 | artifact | computational | Available | 최초 raw-evidence 아티팩트 |
| R5662 | attempt | executable | Partial | Zenodo v2(10.5281/zenodo.21988177) 정정판 |
| R5663 | claim | executable | Supported | 외부 검토 전 계산적 해결 주장 |
| R5664 | artifact | executable | Available | SHA-256 고정 v2 아티팩트 |

이전 작업 지시문에 적힌 `R5747`, `R5748`, `R5749`는 P2624 기록이 아니다. 실시간 페이지 확인 결과 R5747–R5748은 P4468(순환그래프 C(331;1,2)의 독립수), R5749는 P40(Graceful tree conjecture)에 속한다. 따라서 B2 정리 작업에서 이 세 ID를 보존·수정 대상으로 사용하면 안 된다.

## 제3자 독립 감사 R6088–R6090

TheoremDB 기록에 적힌 검증 범위는 다음과 같다.

- 패키지명: `P2624_certified_resolution.zip`
- 패키지 SHA-256: `facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35`
- 13-set 직접 검증: 최대 ordered-difference multiplicity 2
- mod-20 quotient occupancy type: 7개
- admissible raw occupancy vector: 204,360개
- affine orbit: 1,341개
- proof-tree node: 299,903,736개
- 독립 checker: 업로드된 소스에서 C++17로 새로 컴파일
- 음성대조: 첫 proof byte의 예약 비트를 변조한 파일을 거부
- 기록된 결론: 14-set 부존재, 따라서 정확한 최댓값 13

하지만 R6090의 structured packet에는 `source.url = null`, `source.locator = null`로 기록되어 있다. 패키지명과 해시는 공개되어 있으나 원 ZIP의 공개 다운로드 위치는 기록되지 않았다. 2026-08-20 웹 검색에서도 패키지명·SHA·node count에 대응하는 GitHub/Zenodo 공개 원본을 찾지 못했다. 따라서 현재 공개 페이지 정보만으로 제3자가 새로 다운로드하여 재생할 수 있는 불변 아카이브는 확인되지 않는다.

## 정리 원칙

1. 기존 Zenodo v2와 R5659–R5664는 삭제하지 않고 provenance로 유지한다.
2. 새 공개본은 R6088–R6090이 감사한 정확한 바이트의 인증 패키지를 확보하는 경우 그 SHA를 그대로 보존한다.
3. 그 패키지를 확보할 수 없으면 동등한 정적 인증 패키지를 독립 재구성하고 새 SHA·새 검증 로그로 구분한다.
4. 공개 다운로드 URL, 전체 manifest, checker, 생성기 또는 증명 생성 절차, 실행 명령, 음성대조를 함께 제공한다.
5. 그 뒤 TheoremDB에 packet proposal을 제출하고, 검토 전에는 공식 상태를 `Open`에서 임의로 바꾸지 않는다.

## 증거 경계

이 문서는 2026-08-20에 로그인된 TheoremDB 웹 화면과 각 공개 record의 structured packet을 직접 읽어 정리한 것이다. R6088–R6090의 계산 자체를 이 작업공간에서 재실행한 것은 아니며, 공개되지 않은 `P2624_certified_resolution.zip`의 내용도 확인하지 못했다.
