# vHAlign4

[HAlign4](https://academic.oup.com/bioinformatics/article/40/12/btae718/7912339?login=false) is a high-performance multiple sequence alignment software based on the star alignment strategy, designed for efficiently aligning large numbers of sequences. Compared to its predecessor HAlign3, HAlign4 further enhances the ability to handle long sequences and large-scale datasets, enabling fast and efficient alignment on standard computing devices.

## Background
[HAlign3](https://github.com/malabz/HAlign-3) was implemented in Java and was capable of efficiently aligning ultra-large sets of similar DNA/RNA sequences, but had limitations when dealing with long sequences and very large datasets. To address these issues, HAlign4 was reimplemented in C++ and incorporates the Burrows-Wheeler Transform (BWT) and wavefront alignment algorithm.

### Key Improvements
- **Algorithm Optimization**: Replaced the original suffix tree with BWT for more efficient indexing and searching.
- **Memory and Speed Optimization**: Introduced the wavefront alignment algorithm to reduce memory usage and improve alignment speed, especially for long sequences.


## Installation Methods

### 1. Conda Installation

This method uses Conda to create an isolated environment and install HAlign4 with all dependencies.
   
```bash
conda create -n halign4
conda activate halign4
conda install -c malab halign4
halign4 -h
```


### 2. Source Installation with Make

This method allows you to manually compile HAlign4 using Make. For [here](https://github.com/metaphysicser/HAlign-4/releases) to get lastest release.

#### Steps:

1. **Install required packages**:

   You need to have `Make` and a C++ compiler installed. On Ubuntu, you can install them as follows:

   ```bash
   sudo apt-get update
   sudo apt-get install gcc
   ```

2. **Clone the repository**:

   ```bash
   git clone --recursive https://github.com/metaphysicser/HAlign-4.git
   cd HAlign-4
   ```

3. **Compile with Make**:

   ```bash
   make
   ```

4. **Run HAlign4**:

   After the build process completes, the executable `halign4` will be available in the current directory. You can run it using the following command:

   ```bash
   ./halign4 -h
   ```
## Usage
By following either of these methods, you can successfully install and run HAlign4. Choose the method that best fits your development environment.
```
Usage:
  ./vhalign4 <dataset_fasta> <reference_fasta> <output_prefix> [options]

Positional arguments:
  dataset_fasta        Input dataset FASTA file containing query sequences.
  reference_fasta      Reference FASTA file (should contain exactly one sequence).
  output_prefix        Prefix for output files (e.g., output_prefix.fasta / output_prefix.vcf).

Optional arguments:
  -t, --threads <int>  Number of threads to use (default: 1).
  -s, --save-vcf       Enable VCF output generation.
  -h, --help           Show this help message.
```

### Parameter Description
- `Input_file`: Path to the input file or folder (please use `.fasta` as the file suffix or input folder).
- `Output_file`: Path to the output file (please use `.fasta` as the file suffix).
- `-r/--reference`: Reference sequence name (please remove all whitespace), default is the longest sequence.
- `-t/--threads`: Number of threads to use, default is 1.
- `-sa/--sa`: Global `sa` threshold, default is 15.
- `-h/--help`: Show help information.

## Example
Here is a simple example of using HAlign4 for multiple sequence alignment:

```bash
./halign4 input.fasta output.fasta -t 4
```

This command will use 4 threads to align the sequences in the `input.fasta` file, and the result will be saved in the `output.fasta` file.

## License
HAlign4 is developed by the Malab team under the [MIT License](https://github.com/metaphysicser/HAlign-4/blob/main/LICENSE).


