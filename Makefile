.PHONY: build run clean

DOCKER_IMAGE ?= bailey-hp

# Build Docker image that contains compiled binaries under /work/build
build:
	docker build -t $(DOCKER_IMAGE) .

# capture extra words after 'run' (use 'make run -- cmd args...')
RUN_TOKENS := $(filter-out run --,$(MAKECMDGOALS))
RUN_EXE := $(if $(RUN_TOKENS),$(firstword $(RUN_TOKENS)),cg_solver)
RUN_ARGS := $(wordlist 2,$(words $(RUN_TOKENS)),$(RUN_TOKENS))

# Run a built binary inside the Docker image
run: build
	@if [ -n "$(CMD)" ]; then \
	  docker run --rm $(DOCKER_IMAGE) sh -lc "exec /work/build/$(CMD)"; \
	else \
	  docker run --rm $(DOCKER_IMAGE) /work/build/$(RUN_EXE) $(RUN_ARGS); \
	fi

# Clean local CMake build directory (optional)
clean:
	rm -rf build

# swallow extra words after 'run' so make doesn't try to build them as targets
%:
	@:
