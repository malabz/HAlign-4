# vHAlign4

vHAlign4 is a multi-threaded, reference-guided sequence alignment tool built on WFA2 (Wavefront Alignment). It is designed for efficient alignment and downstream processing of viral genomic sequences.

This project is built using CMake and includes WFA2-lib as a submodule.

---

## 1. Build Requirements

### Basic Dependencies

- CMake ≥ 3.22
- C++ compiler with C++17 support
- GCC ≥ 9
- Git (required for fetching submodules)

### Operating System

- Linux

---

## 2. Getting the Source Code (Including Submodules)

⚠️ Important: This project depends on the WFA2-lib submodule, which must be cloned together with the main repository.

```bash
git clone --recursive https://github.com/malabz/HAlign-4.git
```

After cloning, make sure the following directory exists:

```text
PairwiseAlignment/WFA2-lib/
```

---

## 3. Build Instructions (Standard CMake Workflow)

### 1. Create the build directory

```bash
cd HAlign4
mkdir build
cd build
```

### 2. Run CMake configuration

```bash
cmake ..
```

### 3. Compile

```bash
cmake --build . -j
```

Or equivalently:

```bash
make -j
```

---

## 4. Running the Program

After a successful build, the executable will be generated:

```bash
vhalign4
```

Example usage:

```bash
./vhalign4 --help
```

(For detailed options, refer to the program output or subsequent documentation.)

---

## Citation and Acknowledgements
This project makes use of WFA2 (Wavefront Alignment Algorithm):

> Marco-Sola, S., et al. *Fast gap-affine pairwise alignment using the wavefront algorithm*. Bioinformatics.



