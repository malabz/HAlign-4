# HAlign-4

For detailed usage and real-world examples (including the bundled datasets), see: **[`docs/usage.md`](docs/usage.md)**.

 ## Requirements
 ### Build tools

- CMake >= 3.28
- A C++20 compiler (GCC/Clang/MSVC)

### Optional dependencies

- OpenMP (currently REQUIRED by the CMake configuration)
- libcurl (optional: used for downloading URL inputs; if missing, runtime may fall back to shell tools like `curl`/`wget` depending on your environment)
- zlib (optional: used for reading `.gz` FASTA)

---

## Build

This project uses CMake. A typical out-of-source build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The produced binary is:

- `build/halign4`

> Tip
>
> The root `CMakeLists.txt` defaults `CMAKE_BUILD_TYPE` to `Release` if you don’t set it.

---

## Install (Conda)

If you prefer installing via conda (recommended for users who dont want to build from source every time), you can install HAlign-4 from a conda package.

> Note
>
> The exact channel name may differ depending on how you publish the recipe in `recipe/meta.yaml`.
> If you already have a channel (e.g. `bioconda`/`conda-forge`/your own channel), replace the placeholder channel accordingly.

```bash
conda install -c <your-channel> halign4
```

After installation, you should have the `halign4` executable on PATH:

```bash
halign4 --version
```

## Run

Quick start examples and detailed parameter explanations are in **[`docs/usage.md`](docs/usage.md)**.

`halign4` requires two mandatory arguments:

- `-i, --input` : input FASTA path (**must exist**; validated by CLI11 and runtime checks)
- `-o, --output`: output FASTA path

`-w, --workdir` is optional:

- If not provided, a default workdir will be generated: `./tmp-<random>`

Example (with the bundled example data):

```bash
./build/halign4 \
  -i example/data/mt1x.fasta.gz \
  -o out.fasta
```

### Important workdir behavior

In the current implementation, whether a non-empty `workdir` is allowed depends on the build mode (see `checkOption` in `src/halign4.cpp`):

- **Release**: `workdir` must be empty (to avoid overwriting prior results)
- **Debug**: `workdir` may be non-empty (more convenient for iterative development)

### Important CLI options

Defined in `include/config.hpp`:

- `-v, --version`: print version and exit (also shown in `-h/--help`)
- `-t, --thread <int>`: number of threads (default: hardware concurrency)
- `--cons-n <int>`: number of sequences selected in preprocess Top-K (by length) (default: 1000)
- `--kmer-size <int>`: minimizer k-mer size (default: 15; range: 4..31)
- `--kmer-window <int>`: minimizer window size (default: 10)
- `--sketch-size <int>`: Mash/MinHash sketch size (default: 2000)
- `-c, --center-path <path>`: optional; a center/reference FASTA file path (must exist)
- `-p, --msa-cmd <string>`: MSA method keyword or a custom command template string
- `--keep-first-length`: keep only the first/center sequence length unchanged
- `--keep-all-length`: keep all center sequences’ lengths unchanged

### MSA command template (`-p/--msa-cmd`)

#### 1) Default template

The code ships with a built-in default template (see `DEFAULT_MSA_CMD`/`MINIPOA_CMD` in `include/config.hpp`):

```text
minipoa {input} -S -t {thread} -r1 > {output}
```

#### 2) Built-in keywords

Instead of providing a full template, you can pass one of the built-in keywords:

- `-p minipoa` 
 uses the built-in minipoa template
- `-p mafft` 
 uses the built-in MAFFT template
- `-p clustalo` 
 uses the built-in Clustal Omega template

These keywords are expanded internally to templates containing `{input}`, `{output}`, and `{thread}`.

#### 3) Custom template string

You can pass a custom **command template string** via `-p/--msa-cmd`.

Template substitution rules:

- `{input}` and `{output}` are required
- `{thread}` (if present) will be replaced by the thread count

> Note
>
> The program will run a template self-check during argument validation (`cmd::testCommandTemplate`) using a tiny FASTA.
> If your environment doesnt have the selected MSA tool (`minipoa`/`mafft`/`clustalo`) on PATH,
> the self-check will fail early with an error.

---

## Tests

Unit tests are implemented with **doctest** and built via `test/CMakeLists.txt`.

### Recommended: use the helper script

```bash
cd test
./run_tests.sh -t Release -j 8
```

The script will:

- configure and build `build-test/`
- run CTest
- when `--perf` is provided, export `HALIGN4_RUN_PERF=1` (to enable longer perf tests)

Example:

```bash
cd test
./run_tests.sh --perf -- -R align
```

### Run CTest directly

```bash
ctest --test-dir build-test -V
```

Registered test suites currently include (see `test/CMakeLists.txt`):

- `consensus`
- `read_fasta`
- `minimizer`
- `mash`
- `jaccard`
- `write_fasta`
- `align`
- `anchor`

---

## Citation

If you use HAlign-4 in academic work, please cite:

https://doi.org/10.1093/bioinformatics/btae718
