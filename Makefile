.PHONY: build run clean dev-build dev-up dev-down dev-logs dev-exec dev-authorize-key

DOCKER_IMAGE ?= bailey-hp
# Host project root (quoted for safe volume mounts)
HOST_PWD := $(shell pwd)
# Persist outputs generated inside the container
# Also mount inputs read-only for experiments
DOCKER_RUN_FLAGS := -v "$(HOST_PWD)/outputs:/work/outputs" -v "$(HOST_PWD)/inputs:/work/inputs:ro"

# Skip image build if NO_BUILD=1 (run-fast)
NO_BUILD ?= 0

# Build Docker image that contains compiled binaries under /work/build
build:
	docker build -t $(DOCKER_IMAGE) .

# capture extra words after 'run' (use 'make run -- cmd args...')
RUN_TOKENS := $(filter-out run --,$(MAKECMDGOALS))
RUN_EXE := $(if $(RUN_TOKENS),$(firstword $(RUN_TOKENS)),cg_solver)
RUN_ARGS := $(wordlist 2,$(words $(RUN_TOKENS)),$(RUN_TOKENS))

# Run a built binary inside the Docker image
# Always define the rule; conditionally add prerequisite
RUN_DEPS :=
ifneq ($(NO_BUILD),1)
RUN_DEPS += build
endif

run: $(RUN_DEPS)
	@mkdir -p outputs
	@if [ "$(NO_BUILD)" = "1" ]; then \
	  if ! docker image inspect $(DOCKER_IMAGE) >/dev/null 2>&1; then \
	    echo "[ERROR] Docker image '$(DOCKER_IMAGE)' not found. Run 'make build' first."; \
	    exit 1; \
	  fi; \
	fi
	@if [ -n "$(CMD)" ]; then \
	  docker run --rm $(DOCKER_RUN_FLAGS) $(DOCKER_IMAGE) sh -lc "exec /work/build/$(CMD)"; \
	else \
	  docker run --rm $(DOCKER_RUN_FLAGS) $(DOCKER_IMAGE) /work/build/$(RUN_EXE) $(RUN_ARGS); \
	fi

# Clean local CMake build directory (optional)
clean:
	rm -rf build

# swallow extra words after 'run' so make doesn't try to build them as targets
%:
	@:

# -----------------------------
# Dev container via Compose
# -----------------------------
dev-build:
	docker compose build dev

dev-up:
	docker compose up -d --build dev

dev-down:
	docker compose down

dev-logs:
	docker compose logs -f dev

dev-exec:
	docker compose exec dev bash

# ------------------------------------------------------------------
# One-time: add your SSH public key to the dev container (persistent)
# Usage:
#   make dev-authorize-key                 # uses ~/.ssh/id_ed25519.pub
#   make dev-authorize-key PUBKEY=~/.ssh/id_rsa.pub
# ------------------------------------------------------------------
PUBKEY ?= ~/.ssh/id_ed25519.pub
dev-authorize-key:
	@if [ ! -f "$(PUBKEY)" ]; then \
	  echo "[ERROR] Public key not found: $(PUBKEY)"; \
	  echo "Hint: set PUBKEY=~/.ssh/id_rsa.pub or another .pub file"; \
	  exit 1; \
	fi
	@echo "Authorizing key: $(PUBKEY) -> dev@container";
	docker compose exec -T dev bash -lc 'install -d -m 700 -o dev -g dev ~dev/.ssh && cat >> ~dev/.ssh/authorized_keys && chown dev:dev ~dev/.ssh/authorized_keys && chmod 600 ~dev/.ssh/authorized_keys' < $(PUBKEY)
