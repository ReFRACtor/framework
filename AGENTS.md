# AGENTS.md — Coding Agent Guidelines for ReFRACtor Framework

## Project Overview

ReFRACtor (Reusable Framework for Retrieval of Atmospheric Composition) is a multi-instrument Earth science atmospheric composition radiative transfer and retrieval framework. It transforms radiance data from co-located instruments into physical quantities (e.g., ozone VMR) using optimal estimation-based retrieval.

**Core Technologies**: C++14, Fortran 90, Python (via SWIG), CMake

## Build System

### Initial Setup

```bash
# Configure environment and install Python dependencies
./configure.sh

# Install package in editable mode and set up inputs
./python.sh
```

### CMake Build

The project uses CMake with out-of-source builds in the `build/` directory:

```bash
# Configure (done by configure.sh/python.sh typically)
mkdir -p build/debug
cd build/debug
cmake ../.. -DCMAKE_BUILD_TYPE=Debug

# Build
make -j$(nproc)

# Install to build/debug/install by default
make install
```

**Build Types**:
- `Release` (default): Optimized build
- `Debug`: With range checking (`-DBZ_DEBUG` for Blitz++)
- `ReleaseFast`: Extra optimizations (`-Ofast -march=native`)

**Key CMake Options**:
- `-DABSCO_DIR=<path>`: Path to absorption coefficient tables (required for full tests)
- `-DBUILD_PYTHON_BINDING=ON/OFF`: Enable/disable Python SWIG bindings (default: ON)
- `-DCMAKE_INSTALL_PREFIX=<path>`: Installation directory (default: `build/debug/install`)

### Python Package Structure

Python code uses namespace packages under `refractor.*`:
- Source: `python/refractor/framework/`
- Tests: `python/test/`
- Installed via: `pip install --editable .` (development) or `pip install --editable build/debug/bindings/python` (SWIG bindings)

## Testing

### C++ Unit Tests

```bash
cd build/debug
make test_all          # Build test binary
./test/test_all        # Run all C++ unit tests
./test/test_all --help # See test options (Boost.Test framework)
```

Individual tests can be run by test name/suite using Boost.Test selectors.

### Python Tests

```bash
# From repository root
pytest                           # Run all Python tests
pytest python/test/              # Run specific directory
pytest python/test/test_foo.py   # Run specific file
pytest -k "test_pattern"         # Run tests matching pattern
pytest -v                        # Verbose output
```

Configuration in `pytest.ini`. Tests expect:
- `REFRACTOR_INPUT_PATH`: Path to `input/` directory
- `abscodir`: Path to ABSCO tables (set by `python.sh` in `.env`)

### CTest Integration

```bash
cd build/debug
ctest                  # Run all registered tests
ctest -V               # Verbose output
ctest -R pattern       # Run tests matching regex
```

## Architecture

### Core Design Principles

1. **Interface-Driven**: Abstract interfaces (`lib/Interface/`) separate contracts from implementations (`lib/Implementation/`). Allows swapping implementations without modifying existing code.

2. **Auto-Derivatives**: Augmented classes automatically compute Jacobians (∂radiance/∂parameters) alongside forward calculations. Eliminates error-prone hand-coded derivatives and is faster than finite differences.

3. **Unit System**: Arrays carry units; automatic conversions prevent unit mismatch errors.

4. **State Mapping**: The `StateVector` and state mapping system handle optimal estimation parameter updates, supporting both direct and composite mappings.

5. **Multi-Step Retrieval**: Chain multiple NLLS solver configurations, each minimizing different species/spectral regions while holding others fixed. A posteriori info from one step becomes a priori for the next.

### Directory Structure

```
lib/
├── Interface/          # Abstract base classes defining contracts
├── Implementation/     # Concrete implementations of interfaces
├── RadiativeTransfer/  # RT-specific code (LIDORT, 2Stream, etc.)
├── StateMapping/       # State vector mapping system
├── Support/            # Utility classes (units, arrays, logger, etc.)
├── Python/             # Python-specific utilities and configs
└── *.i                 # SWIG interface files

python/refractor/       # Python namespace packages
├── framework/          # Core Python framework code
└── test/               # Python unit tests

test/
├── unit/              # C++ unit tests
├── full/              # Full integration tests
└── pytest/            # Additional pytest configurations

bindings/python/       # Generated SWIG Python bindings (build output)
```

### Key Components

**Radiative Transfer**: Core RT solvers (LIDORT for multi-scatter, 2Stream, etc.) behind abstract `RadiativeTransfer` interface. Fortran implementations wrapped with zero-copy C++ layer.

**NLLS Solvers**: Multiple Non-Linear Least Squares implementations (Connor, Levenberg-Marquardt variants) behind `ConnorSolver` interface. Different solvers make different state assumptions.

**Absorber/Aerosol/Ground**: Atmospheric component models with state vector integration. Each implements appropriate interfaces for RT and retrieval.

**Configuration Layer**: Lua-based configuration system (in separate instrument repos) connects components. Python configuration layer also available via SWIG.

### State System

The retrieval state is managed through:
- `StateVector`: Central parameter vector being optimized
- `StateVectorObserver`: Components that depend on state subscribe to updates
- `SubStateVectorArray`: Maps subset of state to component parameters
- `StateMappingComposite`: Applies transformations (e.g., log-scale) during retrieval

State elements push updates when `StateVector` changes; use `observer_ptr` (non-owning) to avoid circular ownership.

## Common Development Tasks

### Adding a New Implementation

1. Define interface in `lib/Interface/` if not exists
2. Implement in `lib/Implementation/`
3. Add SWIG bindings (`.i` file) for Python access
4. Register with CMake in `lib/CMakeLists.txt`
5. Add unit tests in `test/unit/`
6. Update Python configuration if needed

### Debugging Fortran/C++ Interface Issues

- Check memory layout: Fortran is column-major, C++ (Blitz++) can be either
- Verify array bounds match on both sides
- Use `BZ_DEBUG` (Debug build) to catch array access errors
- Check for inadvertent copies vs. references in wrapper layer

### Working with SWIG Bindings

- Interfaces defined in `.i` files (one per header typically)
- Include `fp_shared_ptr.i` for shared pointer support
- Use `%feature("notabstract")` for pure virtual classes used as interfaces
- Rebuild: `make _swig_wrap` then reinstall Python package
- SWIG Python files generated in `build/*/bindings/python/`

### Unit System

Arrays can carry units (`DoubleWithUnit`, `ArrayWithUnit`):
```cpp
DoubleWithUnit pressure(100.0, "hPa");  // Automatically converts
double in_pa = pressure.convert("Pa").value;  // Explicit conversion
```

Always check unit expectations at interface boundaries.

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

Rules:
- For codebase questions, first run `graphify query "<question>"` when graphify-out/graph.json exists. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).

## Notes

- ABSCO tables are large (~GBs) and distributed separately. `link-absco.sh` creates symlinks from `~/refractor/absco` or `~/refractor-dev/absco`.
- Conda environment recommended for dependencies (HDF5, GSL, Boost, etc.)
- The framework is designed to be extended by instrument-specific adaptations in separate repositories.
- Forked from NASA's RtRetrievalFramework (OCO/OCO-2 heritage)
