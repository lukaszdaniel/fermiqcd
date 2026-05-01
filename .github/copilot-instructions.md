# Copilot instructions for FermiQCD

This repository is a native C++ scientific codebase for lattice QCD. The project is built with plain Makefiles, not with CMake, Meson, Bazel, or any package manager. Use the existing `Makefile` and `Makefile.common` logic when suggesting build or compile changes.

## Key repository structure

- `Libraries/` contains the core FermiQCD and MDP library sources and headers.
- `Examples/` contains many standalone example programs, each with its own `main()`.
- `Converters/` contains conversion utilities for gauge configuration formats.
- `progs/`, `tests/`, `SciDac2007/`, `Vegas/`, `Junk/` contain additional programs, tests and experimental code.
- `Documentation/` and `Doxygen/` contain documentation sources.

## Build and compile

- Root build is triggered by `make` from project root.
- The root `Makefile` delegates to `Libraries`, `progs`, `Examples`, `Converters`, `SciDac2007`, `Vegas`, `Junk`, and `tests`.
- The code is compiled with `g++` or `clang++` and uses `-std=gnu++20` in `Makefile.common`.
- The project standard is C++20; do not suggest newer C++ standards unless the user explicitly requests them.
- Typical flags: `-g -O3 -flto -Wall -Wextra -Werror -pedantic`.
- Library artifacts: `libmdp.a` and `libfermiqcd.a`.
- For single precision, the code can be compiled with `-DUSE_SINGLE_PRECISION`.
- Parallel MPI support is optional via `-DPARALLEL` and `mpiCC` in local `Makefile.common` variants.

## Code conventions and guidance

- Preserve scientific and numerical semantics; this is a high-precision physics codebase.
- Avoid introducing external dependencies, new frameworks, or package-management systems.
- Prefer minimal, localized changes unless the user explicitly asks for a large refactor.
- Do not assume tests or build automation beyond the existing Makefiles.
- Keep the existing style and naming conventions.
- Document public functions, methods, and classes using Doxygen-style comments: `@brief`, `@param`, `@return`, `@note`, etc.

## When answering in this repository

- Respond in the same language as the user.
- Use precise technical language and be concise.
- When recommending build commands, reference `Makefile` and `Libraries/Makefile`.
- If the user asks for a code fix, analyze the relevant `.cpp` or `.h` file and maintain scientific correctness.
- Do not add new top-level tooling such as `CMakeLists.txt`, `vcpkg`, or `conan` unless explicitly requested.
