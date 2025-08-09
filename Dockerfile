# ============================================================================
# Rocky Linux 9 – DD / DQ / QX + Eigen3 + C++17 疎行列環境   (arm64 & x86_64)
# ============================================================================
FROM rockylinux:9

# ---------------------------------------------------------------------------
# 0) 基本開発環境（先に distro-sync）
# ---------------------------------------------------------------------------
RUN dnf -y update --allowerasing && \
    dnf -y install dnf-plugins-core epel-release && \
    dnf config-manager --set-enabled crb && \
    dnf -y groupinstall "Development Tools" && \
    dnf -y install gcc-gfortran gcc-c++ cmake ninja-build \
                   wget tar git eigen3-devel && \
    dnf clean all

# ---------------------------------------------------------------------------
# BLAS/LAPACK ライブラリを別途インストール
# ---------------------------------------------------------------------------
RUN dnf -y install openblas-devel

# ---------------------------------------------------------------------------
# matio と matio-cpp ライブラリをソースからビルド
# ---------------------------------------------------------------------------
RUN set -eux; cd /tmp; \
    # 必要な依存関係をインストール
    dnf -y install hdf5-devel zlib-devel; \
    # matioをソースからビルド
    git clone --depth 1 https://github.com/tbeu/matio.git; \
    cd matio; \
    mkdir build && cd build; \
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=/usr/local \
             -DMATIO_WITH_HDF5=ON; \
    make -j$(nproc) install; \
    cd /tmp && rm -rf matio; \
    # matio-cppをGitHubからクローン&ビルド
    git clone --depth 1 https://github.com/ami-iit/matio-cpp.git; \
    cd matio-cpp; \
    mkdir build && cd build; \
    cmake .. -DCMAKE_BUILD_TYPE=Release \
             -DCMAKE_INSTALL_PREFIX=/usr/local \
             -DMATIOCPP_ENABLE_TESTS=OFF; \
    make -j$(nproc) install; \
    cd /tmp && rm -rf matio-cpp; \
    # 動的ライブラリパスを更新
    ldconfig

# ---------------------------------------------------------------------------
# 1) Bailey ライブラリ (DD, DQ, QX) をビルド
# ---------------------------------------------------------------------------
ENV \
    DDFUN_VER=v03 DDFUN_DIR=/opt/ddfun/fortran \
    DQFUN_VER=v03 DQFUN_DIR=/opt/dqfun/fortran \
    QXFUN_VER=v01 QXFUN_DIR=/opt/qxfun/fortran

RUN set -eux; cd /tmp; \
    for p in ddfun dqfun qxfun; do \
      ver=$(eval echo \$${p^^}_VER); \
      wget -q https://www.davidhbailey.com/dhbsoftware/${p}-${ver}.tar.gz; \
      tar -xzf ${p}-${ver}.tar.gz; \
      cd ${p}-${ver}/fortran; \
      script_name="gnu-complib-${p:0:2}.scr"; \
      if [ ! -f "$script_name" ]; then script_name="gnu-complib.scr"; fi; \
      chmod +x "$script_name"; ./"$script_name"; \
      # ---------- ライブラリが無ければ生成 ---------- \
      if ! ls lib${p}*.a >/dev/null 2>&1; then \
          ar rcs lib${p:0:2}fun.a *.o; \
      fi; \
      dest=$(eval echo \$${p^^}_DIR); \
      mkdir -p "${dest}"; \
      cp lib${p:0:2}*.a *.o *.mod "${dest}/"; \
      cd /tmp; rm -rf ${p}-${ver}*; \
    done

ENV LD_LIBRARY_PATH=${DDFUN_DIR}:${DQFUN_DIR}:${QXFUN_DIR}

# ---------------------------------------------------------------------------
# 2) Fortran→C ラッパをビルド
# ---------------------------------------------------------------------------
COPY interfaces/bailey_wrappers/*.f90 /tmp/
RUN gfortran -J${DDFUN_DIR} -c /tmp/ddfun_cwrap.f90 -I${DDFUN_DIR} && \
    gfortran -J${DQFUN_DIR} -c /tmp/dqfun_cwrap.f90 -I${DQFUN_DIR} && \
    gfortran -J${QXFUN_DIR} -c /tmp/qxfun_cwrap.f90 -I${QXFUN_DIR} && \
    ar rcs ${DDFUN_DIR}/libddwrap.a ddfun_cwrap.o && \
    ar rcs ${DQFUN_DIR}/libdqwrap.a dqfun_cwrap.o && \
    ar rcs ${QXFUN_DIR}/libqxwrap.a qxfun_cwrap.o && \
    rm /tmp/*.f90

# ---------------------------------------------------------------------------
# 3) プロジェクトを投入して CMake ビルド
# ---------------------------------------------------------------------------
WORKDIR /work
COPY . /work
RUN cmake -S . -B build -G Ninja \
        -DQXFUN_DIR=${QXFUN_DIR} -DDQFUN_DIR=${DQFUN_DIR} -DDDFUN_DIR=${DDFUN_DIR} \
        -DCMAKE_BUILD_TYPE=Release && \
    cmake --build build --config Release && \
    mkdir -p outputs

CMD ["/bin/bash"]
