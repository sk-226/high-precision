Zed Remote Development (SSH) with Docker Compose
================================================

目的
- Zedの画面はホスト、言語サーバ/ビルド/依存解決はDocker内で実行して、ローカル未整備による疑似エラー（赤波線）を無くします。
- clangd/fortls/cmake-lsをコンテナ側で動かし、CMakeの`compile_commands.json`を確実に生成・参照します。

同梱内容（このリポで既に実装済み）
- Dockerfileの`dev`ステージ（clangd/fortls/cmake-language-server/sshd/非rootユーザ）。
- Composeサービス`dev`（127.0.0.1:2222でSSH、`/workspace/high-precision`にリポをマウント）。
- Make補助ターゲット：`make dev-up`/`dev-down`/`dev-logs`/`dev-exec`。
- 互換用リンク：`/work -> /workspace/high-precision`（既定の`/work/inputs`がそのまま動く）。
- CMakeLists：`CMAKE_EXPORT_COMPILE_COMMANDS`を常時ON。

前提
- Docker（`docker compose` v2）とZed（Remote via SSH対応版）。
- ホストにSSH鍵（例：`~/.ssh/id_ed25519.pub`）。

初回セットアップ（永続化あり）
1) devコンテナを起動（ホスト側で実行）
- `make dev-up`（devイメージをビルドし、バックグラウンド起動）
- 状態：`docker compose ps`、ログ：`docker compose logs -f dev`

2) 公開鍵をdevユーザに登録（ホスト側で実行）
  - このリポの `compose.yaml` は `/etc/ssh`（ホスト鍵）と `/home/dev/.ssh` を**名前付きボリュームで永続化**しています。よって以下の登録は「初回のみ」でOKです。
- コンテナID：`CID=$(docker compose ps -q dev)`
- `.ssh`作成（所有者/権限込み）：
  `docker compose exec dev bash -lc 'install -d -m 700 -o dev -g dev /home/dev/.ssh'`
- 公開鍵を配置（貼り付け不要のコピー方式）：
  `docker cp ~/.ssh/id_ed25519.pub "$CID":/home/dev/.ssh/authorized_keys`
- 所有権/権限：
  `docker compose exec dev bash -lc 'chown dev:dev /home/dev/.ssh/authorized_keys && chmod 600 /home/dev/.ssh/authorized_keys'`

3) ZedからSSH接続
- Zed → Remote → Connect via SSH → `dev@localhost:2222`
- Open Folder：`/workspace/high-precision`
- 注意：Compose操作（up/down等）は“ホスト側”で行い、Zedのリモート端末内では実行しない。

4) CMake構成（compile_commands.json生成）
- 正しい環境変数名で設定：
  - `export QXFUN_DIR=/opt/qxfun/fortran`
  - `export DQFUN_DIR=/opt/dqfun/fortran`
  - `export DDFUN_DIR=/opt/ddfun/fortran`
- クリーン構成：
  - `rm -rf build`
  - `cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release`
- 確認：`ls build/compile_commands.json`
- 任意（clangdが拾いやすい）：`ln -sf build/compile_commands.json .`

日常運用
- 接続：Zed → `dev@localhost:2222` → `/workspace/high-precision`
- ビルド：`cmake --build build -j`
- 実行：`./build/cg_solver --matrix nos5 --precision dq --tol 1e-15`
  - 既定入力は`/work/inputs`。devイメージには`/work`シンボリックリンク済みのため、`--input-dir`省略でOK。
  - 明示したい場合は`--input-dir inputs`でも可。

停止/再開（ホスト側で）
- 一時停止（状態保持）：`docker compose stop`
- 再開：`docker compose up -d dev`
- 破棄（完全停止）：`make dev-down`（ボリュームは保持。完全初期化したい場合は `docker volume rm high-precision_ssh_host_keys high-precision_dev_ssh` も実行）

ワンライナー（初回の鍵登録を省力化）
- 既定の公開鍵（`~/.ssh/id_ed25519.pub`）を登録：
  - `make dev-authorize-key`
- 別の鍵を使う：
  - `make dev-authorize-key PUBKEY=~/.ssh/id_rsa.pub`

UID/GID
- ホストと所有権を合わせたいとき：
  `UID=$(id -u) GID=$(id -g) docker compose build dev && docker compose up -d dev`

トラブルシュート
- SSH接続不可：`docker compose logs -f dev`でsshdエラーの有無、`/home/dev/.ssh`の権限（dir=700、file=600）を確認。
- LSPが反応しない：`build/compile_commands.json`の存在・更新、Zedの「Restart Language Servers」を確認。
- ライブラリが見つからない（例：`/libqxwrap.a`）：環境変数名の誤りが典型。正しくは`QXFUN_DIR/DQFUN_DIR/DDFUN_DIR`を設定してからCMakeを再構成。

補足（Makeとの使い分け）
- `make build`/`make run`は“ビルド済みイメージでの単発実行”（ciステージ）。
- 開発中はdevコンテナ内で`./build/...`を直接実行するのが自然です。
