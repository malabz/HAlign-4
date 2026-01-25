# HAlign-4 Usage (CLI)

This document explains the command-line arguments of `halign4` and provides runnable examples using the datasets under `test/data/`.

> Notes
>
> - `halign4` needs **three mandatory arguments**: `-i/--input`, `-o/--output`, `-w/--workdir`.
> - Some options are validated by CLI11 at parse time (for example `-i` requires an existing file).
> - `-p/--msa-cmd` is currently validated as an **existing file path** (a template file), **not** a raw command string.

---

## Quick start

Minimal run (uses built-in defaults):

```bash
./build/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o out.fasta \
  -w work
```

---

## Arguments reference

### Required

#### `-i, --input <path>`
Input sequence file path.

- Expected format: FASTA (optionally gzipped, depending on your build)
- Validation: must exist (`CLI::ExistingFile`)

#### `-o, --output <path>`
Output file path. The main result is written here.

- Validation: required (but does not need to exist)

#### `-w, --workdir <path>`
Working directory used for intermediate files.

- The program will create sub-directories under this path (for example `data/`, `temp/`, `result/`).
- **Important**: some builds may require `workdir` to be empty to avoid overwriting previous outputs (see project README).

---

### Optional

#### `-t, --thread <int>`
Number of CPU threads.

- Default: hardware concurrency
- Range check: `1 .. 100000`

#### `--kmer-size <int>`
K-mer size used by minimizer / seeding logic.

- Default: `15`
- Range check: `4 .. 31`

#### `--kmer-window <int>`
Minimizer window size **in number of k-mers**.

- Default: `10`

#### `--cons-n <int>`
Number of sequences selected for consensus step (Top-K by sequence length).

- Default: `1000`

#### `--sketch-size <int>`
Sketch size used by Mash/MinHash-related components.

- Default: `2000`

#### `-c, --center-path <path>`
Provide an explicit center/reference sequence file (FASTA).

- If provided, the program will use these sequences as the reference/center set instead of auto-selecting.
- Validation: must exist (`CLI::ExistingFile`)

#### `-p, --msa-cmd <path>`
Path to a **text file** that contains the MSA command template.

- Validation: must exist (`CLI::ExistingFile`)
- The template may contain placeholders:
  - `{input}`: required, replaced with the input file path for MSA
  - `{output}`: required, replaced with the output file path for MSA
  - `{thread}`: optional, replaced with the integer from `-t/--thread`

Built-in default template (used when you don’t pass `-p`):

```text
minipoa {input} -S -t {thread} -r1 > {output}
```

Important note:

- In `halign4`, **minipoa is the default high-quality aligner** (via the built-in template above).
- In the examples below we use **MAFFT** only because it’s a common MSA tool and its CLI is easy to demonstrate.

Security note:

- The command is executed via the system shell. Treat template inputs as trusted data.

#### `--keep-first-length`
Keep the **first reference sequence** (i.e. the first record in `-c/--center-path`) ungapped in the final MSA.

- Source-level meaning (matches `RefAligner` implementation):
  - When set, the pipeline removes alignment columns where the reference has a gap, so the first reference sequence will not contain inserted `-` columns.
- When `-c` contains multiple reference sequences:
  - only the **first** reference sequence is guaranteed to have no inserted gaps.

#### `--keep-all-length`
Keep **all reference sequences** in `-c/--center-path` ungapped in the final MSA.

- Source-level meaning (matches `RefAligner` implementation):
  - When set, the pipeline removes alignment columns that would introduce gaps into any of the reference sequences.
- When `-c` contains multiple reference sequences:
  - **every** reference sequence is guaranteed to have no inserted gaps.

> Important clarification
>
> These two flags are about **reference sequences from `-c`** (the “center/reference FASTA”), not about general query sequences.

#### `--save-workdir`
Keep the working directory after successful completion.

- Default behavior: remove workdir on success

---

## Examples

### Example 1: Minimal dataset (`mt1x`) + demonstrate `-p/--msa-cmd` (using MAFFT)

Dataset:

- `test/data/mt1x.fasta.gz`

Goal:

- Show the most basic run.
- Show how to use `-p` correctly (pass a **template file path**).

1) Create a template file (example: `mafft.template.txt`):

```text
mafft --thread {thread} --auto {input} > {output}
```

Why this looks different from the built-in minipoa template:

- MAFFT writes alignment to stdout by default, so redirection (`> {output}`) is the simplest portable pattern.

2) Run:

```bash
./build/halign4 \
  -i test/data/mt1x.fasta.gz \
  -o mt1x.out.fasta \
  -w mt1x.work \
  -t 8 \
  -p mafft.template.txt
```

If you don’t have `mafft` installed, either install it or switch the template file back to a command that exists in your environment.

---

### Example 2: COVID dataset + demonstrate `-c`, `--keep-first-length`, `--keep-all-length`

Dataset:

- `test/data/covid-ref.fasta.gz`
- `test/data/covid-test.fasta.gz`

Background:

- The first record in `covid-ref.fasta.gz` is the Wuhan reference sequence.
- Other sequences are consensus sequences for variants produced by:
  https://github.com/corneliusroemer/pango-sequences

#### 2.1 Keep only the *first* reference sequence ungapped (`--keep-first-length`)

```bash
./build/halign4 \
  -i test/data/covid-test.fasta.gz \
  -o covid.keep_first.out.fasta \
  -w covid.keep_first.work \
  -c test/data/covid-ref.fasta.gz \
  --keep-first-length
```

#### 2.2 Keep *all* reference sequences ungapped (`--keep-all-length`)

```bash
./build/halign4 \
  -i test/data/covid-test.fasta.gz \
  -o covid.keep_all.out.fasta \
  -w covid.keep_all.work \
  -c test/data/covid-ref.fasta.gz \
  --keep-all-length
```

---

## How to understand `--keep-first-length` vs `--keep-all-length` (toy example)

This section uses a tiny, **made-up** example to make the idea concrete. All alignment rows below have the **same length**.

Assume your `-c/--center-path` has **two** reference sequences:

```text
ref1 (first):  ACGTAC
ref2:          ACGTC
```

And you have one query sequence:

```text
q1:            ACGTTAC
```

A generic MSA tool may produce an alignment like this (gaps inserted into references are allowed):

```text
ref1  ACGT-AC
ref2  ACGT--C
q1    ACGTTAC
```

### Case A: default behavior (no `--keep-...` flags)

- Gaps in reference sequences are allowed.
- This can change the *effective* coordinate system of references.

### Case B: `--keep-first-length`

- Guarantee: **the first reference sequence (`ref1`) will not contain inserted gaps**.
- But other references (like `ref2`) may still have gaps.

Conceptually, the pipeline will drop columns where `ref1` is `-`, so the alignment becomes:

```text
ref1  ACGTAC   (ref1 has no '-')
ref2  ACGT-C   (column was removed, ref2 becomes ungapped here)
q1    ACGTAC   (query is projected accordingly)
```

### Case C: `--keep-all-length`

- Guarantee: **all references in `-c` will not contain inserted gaps**.
- So the pipeline will drop columns that would introduce gaps in any reference.

Using the original MSA snippet above, dropping the column where `ref1` has `-` *and* the column where `ref2` has `-` yields:

```text
ref1  ACGTC
ref2  ACGTC
q1    ACGTC
```

What to take away:

- Both flags can reduce the number of alignment columns (because they remove “reference-gap columns”).
- `--keep-all-length` is stricter, so it may remove more columns than `--keep-first-length` when `-c` contains multiple reference sequences.

> This toy example is meant to build intuition; exact output depends on the real sequences and the chosen MSA method.

---

## Troubleshooting

- **`-p/--msa-cmd` rejects my command string**: this is expected. `-p` needs a *file path* to a template file (validated by `ExistingFile`).
- **Workdir already exists**: remove it, choose a new `-w`, or build/run in Debug mode if your build allows reusing a non-empty workdir.
