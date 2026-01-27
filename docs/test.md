# Tests

This document explains how to run HAlign-4 tests under the `test/` directory.

HAlign-4 uses **doctest** for unit tests and CTest for registration.

---

## 1. Quick start (recommended)

Use the helper script:

```bash
cd test
./run_tests.sh -t Release -j 8
```

What it does:

- configures the test build into `build-test/` (or your chosen `-B` directory)
- builds the test binaries
- runs `ctest` (verbose by default)

---

## 2. Common workflows

### 2.1 Run a single CTest test by name

Example: run the `align` test only:

```bash
cd test
./run_tests.sh -t Release -- -- -R align
```

(Everything after `--` is passed to `ctest`.)

### 2.2 Enable perf tests

Some tests are heavier and are guarded by an environment variable.

```bash
cd test
./run_tests.sh -t Release --perf
```

This sets:

- `HALIGN4_RUN_PERF=1`

### 2.3 Build directory management

To build tests into a custom directory:

```bash
cd test
./run_tests.sh -B ../build-test-release -t Release -j 16
```

To clean rebuild:

```bash
cd test
./run_tests.sh --clean -t Release
```

---

## 3. Run CTest directly (advanced)

If you already have a configured build directory:

```bash
ctest --test-dir build-test -V
```

---

## 4. Using `test/data` as example inputs

The repository includes small datasets that are useful for smoke tests and documentation examples:

- `test/data/mt1x.fasta.gz` (minimal dataset)
- `test/data/covid-ref.fasta.gz` / `test/data/covid-test.fasta.gz`

For runnable CLI examples, see [`docs/usage.md`](usage.md).

