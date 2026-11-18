#!/usr/bin/env bash
# Fetch and build a serial shared liblammps into .lammps-install, keeping the
# source headers (include-src = src/) that the direct-internals embed
# compiles against. This file is the only copy of the release tag. CI reads
# it with --print-tag, exports LAMMPS_TAG, and uses that same value as the
# cache key.
set -euo pipefail

# Source of truth for the liblammps release.
DEFAULT_LAMMPS_TAG=stable_22Jul2025_update2
LAMMPS_TAG="${LAMMPS_TAG:-${DEFAULT_LAMMPS_TAG}}"

if [ "${1:-}" = "--print-tag" ]; then
  printf '%s\n' "${LAMMPS_TAG}"
  exit 0
fi

PREFIX="${PWD}/.lammps-install"

if [ -f "${PREFIX}/.built-${LAMMPS_TAG}" ]; then
  echo "liblammps ${LAMMPS_TAG} already built at ${PREFIX}"
  exit 0
fi

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir:?}"' EXIT

curl -sL "https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_TAG}.tar.gz" \
  -o "${workdir}/lammps.tar.gz"
tar -C "${workdir}" -xzf "${workdir}/lammps.tar.gz"
src_dir="${workdir}/lammps-${LAMMPS_TAG}"

cmake -S "${src_dir}/cmake" -B "${workdir}/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MPI=no \
  -DBUILD_OMP=no \
  -DPKG_MANYBODY=on \
  -DBUILD_SHARED_LIBS=on \
  -DBUILD_TOOLS=no \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}"
cmake --build "${workdir}/build" -j "$(nproc)"
cmake --install "${workdir}/build"

mkdir -p "${PREFIX}/include-src"
cp "${src_dir}/src/"*.h "${PREFIX}/include-src/"
mkdir -p "${PREFIX}/include-src/STUBS"
cp "${src_dir}/src/STUBS/"*.h "${PREFIX}/include-src/STUBS/"
mkdir -p "${PREFIX}/include-src/fmt"
cp "${src_dir}/src/fmt/"*.h "${PREFIX}/include-src/fmt/"
mkdir -p "${PREFIX}/include-src/nlohmann"
cp "${src_dir}/src/nlohmann/"*.hpp "${PREFIX}/include-src/nlohmann/"

touch "${PREFIX}/.built-${LAMMPS_TAG}"
echo "liblammps ${LAMMPS_TAG} installed at ${PREFIX}"
