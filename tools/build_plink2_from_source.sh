#!/usr/bin/env bash
# Documents the exact way plink2 was built for the V5 run reported in TEST_RESULTS.md.
# Only needed if you cannot download an official binary. Prefer the official release binary.
# Build environment had no BLAS/LAPACK headers (cblas.h missing), hence NO_LAPACK=1.
set -euo pipefail
git clone https://github.com/chrchang/plink-ng.git
cd plink-ng
git checkout a2b85291f90a1c6090bd9039d47c3314176768b2   # commit used; reported as "PLINK v2.0.0-b.1-devNL"
cd 2.0/build_dynamic
# With NO_LAPACK=1 this commit failed to compile on a missing fabs(); the one-line patch below
# adds <math.h>. If your checkout already compiles, skip it.
grep -q '#include <math.h>' ../plink2_matrix.h 2>/dev/null || \
  sed -i '0,/#include/s//#include <math.h>\n#include/' ../plink2_matrix.h
make -j"$(nproc)" NO_LAPACK=1
echo "Built: $(pwd)/plink2   ->   export PLINK2=$(pwd)/plink2"
# NOTE: NO_LAPACK=1 disables LAPACK-backed features (including plink2 --pca). Script 02's
# default python PCA engine does not need it; --pca-engine plink2 does and is UNTESTED here.
