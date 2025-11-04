## TODO: 高精度型の出力精度制御と数値特性の明確化（DD/DQ/QX）

本ドキュメントは、Bailey 系高精度型（`bailey::DDNumber`/`DQNumber`/`QXNumber`）の
標準入出力および数値特性の扱いを改善し、実装時の手戻りを最小化するための詳細設計です。

### 背景と目的
- 現状の `operator<<` は固定桁で `*_to_string` を呼んでおり、`std::cout << std::setprecision(N)` の設定を反映しません。
- `std::numeric_limits<T>::max_digits10` を持たないため、型ごとの「ラウンドトリップ安全な出力桁数」を簡単に参照できません。
- 目的は次の2点:
  - 出力時にストリーム精度（`os.precision()`）を尊重し、適切に桁数を制御する。
  - `std::numeric_limits<T>` をユーザー定義型に特殊化し、`digits`/`digits10`/`max_digits10`/`epsilon` などの特性を明示する。

### スコープ
- 対象: `bailey::DDNumber`, `bailey::DQNumber`, `bailey::QXNumber`
- 非対象: アルゴリズム本体（CG 等）の数値改善。今回は I/O と数値特性の宣言に限定。

### 設計方針

#### 1) `std::numeric_limits<T>` の特殊化を用意
- 目的: `max_digits10` を中心とした数値特性を標準 API から取得できるようにする。
- 配置案: `include/bailey/numeric_limits.hpp` を新設し、3型の完全特殊化を定義。
- 必須フィールド（最低限）:
  - `static constexpr bool is_specialized = true;`
  - `static constexpr int radix = 2;`
  - `static constexpr int digits;`（2進仮数部ビット数）
  - `static constexpr int digits10;`（10進有効桁）
  - `static constexpr int max_digits10;`（ラウンドトリップ用安全桁）
  - `static constexpr T epsilon() noexcept;`（≒ `2^(1-digits)`）
- 推奨値（理論/既存宣言と整合）:
  - DD: `digits=106`, `digits10=31`, `max_digits10=33`, `epsilon=ldexp(1.0L, -105)`
  - DQ: `digits=212`, `digits10=66`, `max_digits10=65`, `epsilon=ldexp(1.0L, -211)`
  - QX: `digits=113`, `digits10=33`, `max_digits10=36`, `epsilon=ldexp(1.0L, -112)`
    - 備考: `max_digits10 = ceil(1 + digits*log10(2))`。53→17, 113→36 に一致。
- 注意点:
  - `std` 名前空間へのユーザー定義型の特殊化は標準で許可されています。
  - `highest()/lowest()/infinity()` 等は必要になった時点で拡張（当面は未使用）。

#### 2) 出力演算子 `operator<<` でストリーム精度を尊重
- 現状: 固定の `digits` を `*_to_string` に与えている。
- 変更: `os.precision()` を読み、その値を `1..numeric_limits<T>::max_digits10` でクランプして渡す。
- 実装ポイント:
  - バッファ長を余裕あるサイズに（例: DD: 80→128、DQ/QX: 128 以上）。
  - `#include <algorithm>` を入れて `std::clamp` を利用。
  - フォールバック: もし `numeric_limits<T>::max_digits10` が未定義なビルド環境では、型ごとの安全既定値を使う。

スニペット（DD の例、DQ/QX も同様に）：
```cpp
// dd_arithmetic.hpp
#include <algorithm>
#include <limits>

inline std::ostream& operator<<(std::ostream& os, const bailey::DDNumber& d) {
    char s[128] = {0};
    int req = static_cast<int>(os.precision());
    int maxd = std::numeric_limits<bailey::DDNumber>::max_digits10; // 33
    int digits = std::clamp(req, 1, maxd);
    dd_to_string(d.dd, &digits, s, sizeof(s));
    return os << s;
}
```

#### 3) `Eigen::NumTraits<T>` の微修正（任意だが推奨）
- `using Literal = ...;` を明示しておくと、Eigen 内で生成される数値リテラルの型が明確になります。
  - DD/DQ: `using Literal = double;`
  - QX: `using Literal = long double;`
- 既存の `digits()`/`digits10()` は維持。

### 既存コードへの影響
- API 互換: 破壊的変更なし。`operator<<` の出力桁数が「固定」→「ストリーム設定尊重」に変わるのみ。
- ログ差分: 出力の既定桁は実行時の `std::cout.precision()` に依存（既存テストは `std::setprecision` を明示しているため影響軽微）。

### 実装手順（作業順）
1. `include/bailey/numeric_limits.hpp` を新規作成し、3型の `std::numeric_limits` 特殊化を定義。
2. `include/bailey/dd_arithmetic.hpp`/`dq_arithmetic.hpp`/`qx_arithmetic.hpp` の `operator<<` を、`os.precision()` を尊重する形に変更。
   - 必要なら一時的に `digits` 既定値を `max_digits10` に設定。
3. `Eigen::NumTraits` に `using Literal = ...;` を追加（任意）。
4. ビルドして `tests/eigen_precision_test.cpp` を実行。
   - `std::setprecision(std::numeric_limits<Scalar>::max_digits10)` で各型の表示が適正桁で出ることを確認。

### テスト計画
- 単体テスト（新規）: `tests/precision_io_test.cpp`（仮）
  - 各型について、`std::setprecision(std::numeric_limits<T>::max_digits10)` を設定し、`to_string` 長または出力に含まれる有効数字の概算が期待値±1以内であること。
  - 代表値（0, 1, π, √2, 極小/極大近辺）を対象。
- 既存 `tests/eigen_precision_test.cpp` の最小変更: `std::setprecision(max_digits10)` を使用するサンプルを追記（比較はしない）。

### リスクと対策
- `std::numeric_limits` の特殊化は ODR に注意: 1 つのヘッダにまとめ、重複定義を避ける。
- Fortran 側 `*_to_string` の桁数上限: ライブラリが内部でクランプする場合があるため、`digits <= max_digits10` を厳守。
- バッファ長不足: 余裕あるサイズ（>= 2 + max_digits10 + 余白）を確保。

### 受け入れ条件（Acceptance Criteria）
- `std::cout << std::setprecision(std::numeric_limits<T>::max_digits10)` を使うと、各型で希望桁数の文字列が出力される。
- 既存の演算・テストは非変更でビルド・実行が通る。
- `digits`/`digits10`/`max_digits10`/`epsilon` が `std::numeric_limits<T>` で取得できる。

### 参考実装の雛形

`include/bailey/numeric_limits.hpp`（抜粋）:
```cpp
#pragma once
#include <limits>
#include <cmath>
#include "bailey/dd_arithmetic.hpp"
#include "bailey/dq_arithmetic.hpp"
#include "bailey/qx_arithmetic.hpp"

namespace std {
template<> struct numeric_limits<bailey::DDNumber> {
    static constexpr bool is_specialized = true;
    static constexpr int radix = 2;
    static constexpr int digits = 106;
    static constexpr int digits10 = 31;
    static constexpr int max_digits10 = 33;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bailey::DDNumber epsilon() noexcept {
        return bailey::DDNumber(static_cast<double>(ldexp(1.0, -105)));
    }
};

template<> struct numeric_limits<bailey::DQNumber> {
    static constexpr bool is_specialized = true;
    static constexpr int radix = 2;
    static constexpr int digits = 212;
    static constexpr int digits10 = 66;
    static constexpr int max_digits10 = 65;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bailey::DQNumber epsilon() noexcept {
        long double e = std::ldexp(1.0L, -211);
        return bailey::DQNumber(static_cast<double>(e));
    }
};

template<> struct numeric_limits<bailey::QXNumber> {
    static constexpr bool is_specialized = true;
    static constexpr int radix = 2;
    static constexpr int digits = 113;
    static constexpr int digits10 = 33;
    static constexpr int max_digits10 = 36;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bailey::QXNumber epsilon() noexcept {
        return bailey::QXNumber(std::ldexp(1.0L, -112));
    }
};
} // namespace std
```

`operator<<` の変更（DD の例）:
```cpp
inline std::ostream& operator<<(std::ostream& os, const bailey::DDNumber& d) {
    char s[128] = {0};
    int req = static_cast<int>(os.precision());
    int maxd = std::numeric_limits<bailey::DDNumber>::max_digits10;
    int digits = std::clamp(req, 1, maxd);
    dd_to_string(d.dd, &digits, s, sizeof(s));
    return os << s;
}
```

### フォローアップ（任意強化）
- `Eigen::NumTraits<T>` に `using Literal = ...;` を追加し、定数生成の一貫性を向上。
- `std::numeric_limits<T>` に `min_exponent`/`max_exponent`/`lowest`/`highest` を段階的に追加。
- `tests/scalar_precision_test.cpp` の QX 許容値を `2^-112` に寄せる（別件）。

### 作業チェックリスト
- [ ] `include/bailey/numeric_limits.hpp` を作成して 3 型の特殊化を実装
- [ ] `dd_arithmetic.hpp`/`dq_arithmetic.hpp`/`qx_arithmetic.hpp` の `operator<<` を修正
- [ ] 必要な `#include <algorithm>`, `#include <limits>` を追加
- [ ] 既存テストをビルド・実行（`eigen_precision_test` で出力桁制御を確認）
- [ ] （任意）`NumTraits` に `Literal` を追加
- [ ] （任意）I/O と特性について `doc/` に簡易ガイドを追加
