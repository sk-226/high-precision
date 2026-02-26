---
name: precision_cast implementation v2
overview: Bailey型間の変換を文字列経由に統一（反復中のベクトル変換も含む）し、DQ→DDなどで double へ落ちる精度損失を解消する。digits は変換元(From)のdecimal_digits()を使用。
todos:
  - id: rev-add-precision-cast
    content: precision_traits.hpp に include追加、from_string<double> と precision_cast（digits=From）を追加
    status: pending
  - id: rev-update-ssor-pcg
    content: ssor_pcg_mixed_precision.hpp の変換を precision_cast に置換（65/75/96/148/157行付近）
    status: pending
  - id: rev-sanity-check
    content: ビルドと簡単な実行で動作確認（性能劣化は注記）
    status: pending
---

# 文字列経由 `precision_cast`（ベクトル変換も含めて常用）

## 目的

`include/algorithms/ssor_pcg_mixed_precision.hpp` の T→P / P→T 変換で `to_double()` に落ちてしまう問題を解消し、例えば **DQ→DD でDDの30桁精度を活かす**。ユーザー要望により **反復中のベクトル変換でも常に文字列経由** を採用する。

## 重要なトレードオフ（明示）

- **精度**: 文字列経由は DQ→DD でも double へ丸めないため精度が保たれる。
- **性能**: 各反復で `n` 要素ぶん `to_string/from_string` を回すため、サイズが大きい問題ではボトルネックになりうる。

## 実装方針

### 1) [include/bailey/precision_traits.hpp](include/bailey/precision_traits.hpp)

#### 追加する include

- `<type_traits>`（`std::is_same_v` 用）
- `<string>`（`std::stod` / `std::string` 生成用。既に他ヘッダ経由で入る可能性はあるが、明示する）

#### 追加するAPI（`namespace bailey`）

- `from_string<double>` の特殊化（`std::stod`）
- `precision_cast<To>(From)`
- `To==From` → そのまま
- `To==double` → `to_double(val)`
- `From==double` → `To(val)`（DD/DQ/QX は `double` コンストラクタあり）
- それ以外（Bailey型同士）→ 文字列経由

**digitsポリシー（ユーザー指定）**

- 文字列化の `digits` は **変換元 `From` の `PrecisionTraits<From>::decimal_digits()`** を使う。
- 例: DQ→DD の場合 `66` 桁で文字列化してから DDに parse（DD 側で丸められる）。

実装イメージ（要点）:

- `auto s = to_string(val, PrecisionTraits<From>::decimal_digits());`
- `return from_string<To>(s);`

追加の安全策:

- 想定外の型が来たときに分かりやすく落とすため `static_assert`（例: double/DD/DQ/QX のみ許可）を入れる。

### 2) [include/algorithms/ssor_pcg_mixed_precision.hpp](include/algorithms/ssor_pcg_mixed_precision.hpp)

`P(to_double(...))` や `T(to_double(...))` を `bailey::precision_cast<>()` に置換。置換対象:

- **65行目** `diag[i] = ...`
- **75行目** `P val = ...`（DL構築）
- **96行目** `P val = ...`（DU構築）
- **148行目** `convert_to_prec` の各要素代入
- **157行目** `convert_from_prec` の各要素代入

## 期待する効果

- `convert_to_prec` が `P(v[i])` のような「未定義な直接キャスト」で落ちる問題を回避。
- DQ→DD で double 丸めを避け、DDの桁を維持。