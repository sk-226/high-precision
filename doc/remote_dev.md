# Remote Development over SSH

This guide describes an editor-agnostic workflow for developing inside the dev container via SSH.

## Prerequisites

- Docker and Docker Compose v2
- OpenSSH client (`ssh`)
- Make (optional, recommended)

## Overview

You will run a persistent dev container exposing SSH on `localhost:2222` and connect as `dev`. The project root inside the container is `/work`.

> [!IMPORTANT]
> Run Docker/Compose commands on the host terminal, not inside the SSH session.

## Start the Dev Container

```bash
make dev-up     # build and start the dev service
```

The service keeps SSH host keys and `authorized_keys` in persistent volumes so you don't need to reconfigure after restarts.

> [!TIP]
> `make dev-logs` tails the service logs; `make dev-down` stops and removes it.

## Authorize Your SSH Key (one-time)

Use the helper target (defaults to `~/.ssh/id_ed25519.pub`):

```bash
make dev-authorize-key
```

Or specify another key:

```bash
make dev-authorize-key PUBKEY=~/.ssh/id_rsa.pub
```

## Connect from Your Editor or Terminal

- SSH target: `dev@localhost` port `2222`
- Project folder to open: `/work`

Examples:

```bash
ssh -p 2222 dev@localhost
# then inside the container:
cd /work
```

## Build and Run Inside the Container

```bash
cd /work

# LSP only: configure project (generates build/compile_commands.json; no binaries)
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug

# Build binaries (optional for LSP)
cmake --build build -j

# Example run
./build/cg_solver --matrix nos7 --precision double --tol 1e-12 --max-iter 30
```

Language servers typically read `build/compile_commands.json`. If your editor expects it at the repo root:

```bash
ln -sf build/compile_commands.json ./compile_commands.json
```

## Inputs and Outputs

- Host `./inputs` is mounted to `/work/inputs` (read-only).
- Host `./outputs` is mounted to `/work/outputs` (read–write) and persists results.

If `./inputs` is empty or missing, matrix loading will fail. Add Matrix Market `.mtx` files under `inputs/`.

## Stop / Restart

```bash
docker compose stop        # pause
docker compose up -d dev   # resume
make dev-down              # remove service (volumes persist)
```
