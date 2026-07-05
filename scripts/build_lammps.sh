#!/usr/bin/env bash
# Fetch and build a serial shared liblammps into .lammps-install, keeping the
# source headers (include-src = src/) that the direct-internals embed
# compiles against. CI caches the whole install dir keyed on LAMMPS_TAG.
set -euo pipefail

LAMMPS_TAG="${LAMMPS_TAG:-stable_22Jul2025}"
PREFIX="${PWD}/.lammps-install"

if [ -f "${PREFIX}/.built-${LAMMPS_TAG}" ]; then
  echo "liblammps ${LAMMPS_TAG} already built at ${PREFIX}"
  exit 0
fi

workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT

curl -sL "https://github.com/lammps/lammps/archive/refs/tags/${LAMMPS_TAG}.tar.gz" \
  -o "${workdir}/lammps.tar.gz"
tar -C "${workdir}" -xzf "${workdir}/lammps.tar.gz"
src_dir="$(find "${workdir}" -maxdepth 1 -type d -name 'lammps-*' | head -1)"

cmake -S "${src_dir}/cmake" -B "${workdir}/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_MPI=no \
  -DBUILD_OMP=no \
  -DBUILD_SHARED_LIBS=on \
  -DBUILD_TOOLS=no \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}"
cmake --build "${workdir}/build" -j "$(nproc)"
cmake --install "${workdir}/build"

mkdir -p "${PREFIX}/include-src"
cp "${src_dir}/src/"*.h "${PREFIX}/include-src/"
mkdir -p "${PREFIX}/include-src/STUBS"
cp "${src_dir}/src/STUBS/"*.h "${PREFIX}/include-src/STUBS/"

touch "${PREFIX}/.built-${LAMMPS_TAG}"
echo "liblammps ${LAMMPS_TAG} installed at ${PREFIX}"
