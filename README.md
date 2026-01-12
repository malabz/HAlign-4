# HAlign-4



## Requirements

### Build tools

- CMake >= 3.28
- A C++20 compiler (GCC/Clang/MSVC)

### Optional dependencies

- OpenMP (required by current CMake)
- libcurl (optional; if missing, downloads fall back to `curl`/`wget` shell tools at runtime)
- zlib (optional; enables gz FASTA input reading and compression support)

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

## Run

`halign4` requires three mandatory arguments:

- `-i, --input` : input FASTA path (currently checked as an existing file by CLI11)
- `-o, --output`: output FASTA path
- `-w, --workdir`: working directory for intermediate files

Example (with the bundled example data):

```bash
./build/halign4 \
  -i example/data/mt1x.fasta.gz \
  -o out.fasta \
  -w work
```

### Important CLI options

Defined in `include/config.hpp`:

- `-t, --thread <int>`: number of threads (default: hardware concurrency)
- `--cons-n <int>`: number of longest sequences selected for the consensus/MSA stage (default: 1000)
- `--kmer-size <int>`: k-mer size for downstream steps (default: 15; constraints 4..31)
- `--kmer-window <int>`: minimizer window size (default: 10)
- `--sketch-size <int>`: sketch size (default: 2000)
- `-c, --center-path <path>`: optional, manually supply a center FASTA record file used for alignment
- `-p, --msa-cmd <path>`: path to an MSA command template file (see below)
- `--keep-first-length`: keep only the first/center sequence length unchanged
- `--keep-all-length`: keep all center sequences’ lengths unchanged

### MSA command template

The default command template is:

```text
minipoa {input} -S -t {thread} -r1 > {output}
```

The tool performs a small self-check by running the template on a tiny FASTA during option validation.

Notes:

- The template **must include** `{input}` and `{output}`.
- If it includes `{thread}`, it will be replaced by the thread count.
- The execution uses the system shell (`std::system`), so treat templates as trusted input.


---

## Tests

Unit tests are implemented with **doctest** and built via `test/CMakeLists.txt`.

You can run tests using the helper script:

```bash
cd test
./run_tests.sh -t Release -j 8
```

Or run CTest directly after configuring `build-test`:

```bash
ctest --test-dir build-test -V
```

Registered test suites include:

- `consensus`
- `read_fasta`
- `minimizer`
- `mash`
- `jaccard`

---



## Citation

If you use HAlign-4 in academic work, please consider adding a citation section once a publication/DOI is available.

