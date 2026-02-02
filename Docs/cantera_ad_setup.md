# Cantera-AD + SU2 integration notes

This branch uses the Cantera-AD fork as a git submodule:
- https://github.com/bigfooted/cantera
- branch: feature_differentiation (pinned by git submodule commit)

## Build Cantera-AD (inside SU2)
From SU2 root:

    cd subprojects/cantera
    rm -rf install_ad cantera.conf
    scons --config=force build -j"$(sysctl -n hw.ncpu)" \
      boost_inc_dir=/opt/homebrew/include \
      prefix="$PWD/install_ad" \
      libdirname="$PWD/install_ad/lib" \
      f90_interface=n \
      system_eigen=n system_fmt=n system_yamlcpp=n \
      system_sundials=n system_blas_lapack=n \
      optimize=n
    scons install prefix="$PWD/install_ad" libdirname="$PWD/install_ad/lib"

## Build SU2 with Cantera enabled
From SU2 root:

    rm -rf build
    ./meson.py setup build -Denable-cantera=true --prefix="$PWD"
    ./ninja -C build
    ./ninja -C build install

Notes:
- `install_ad/` is intentionally not committed.
- Meson subproject wrapper is `subprojects/cantera/meson.build`.
