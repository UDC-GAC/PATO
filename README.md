# <img src=".github/PATO.gif"  height="64pt" style="margin-bottom: -28%"> PATO: high PerformAnce TriplexatOr

PATO: high PerformAnce TriplexatOr is a high performance tool for the fast and
efficient detection of triple helices and triplex features in nucleotide
sequences. PATO: high PerformAnce TriplexatOr is a modern alternative to
[Triplexator](https://github.com/Gurado/triplexator) and
[TDF](https://github.com/CostaLab/reg-gen) and functions as a drop in
replacement to accelerate the triplex analyses in multicore computers. Moreover,
PATO: high PerformAnce TriplexatOr's efficiency allows a more exhaustive search
of the triplex-forming solution space, and so it achieves higher levels of
prediction accuracy in far less time than other tools in the state of the art.

## Compiling PATO

Download the source code from this repository, either use Git or download a copy
from GitHub, and let CMake compile PATO: high PerformAnce TriplexatOr for you:

```bash
$ cmake -B build . && cmake --build build
```

Note that macOS users must explicitly specify an OpenMP-enabled compiler to
compile PATO: high PerformAnce TriplexatOr. For instance, to use `g++-12`
(available via Homebrew), execute:

```bash
$ cmake -B build -D CMAKE_CXX_COMPILER=g++-12 . && cmake --build build
```

## Executing PATO

Now that PATO: high PerformAnce TriplexatOr has been compiled, execute the
application as follows:

```bash
$ ./build/tools/PATO/PATO [options] {-ss tfo_file | -ds tts_file | -ss tfo_file -ds tts_file}
```

Execute `./build/tools/PATO/PATO --help` for a detailed list of execution modes,
options, and flags.

### Setting the number of threads

PATO: high PerformAnce TriplexatOr uses OpenMP to parallelize its triplex search
algorithms. The OpenMP runtime will automatically spawn as many threads as there
are available CPU cores. To reduce the number of threads spawned by the
application one has to explicitly set the `OMP_NUM_THREADS` environment variable
to a value greater than 0. For instance, to run PATO: high PerformAnce
TriplexatOr with 4 threads, execute:

```bash
$ OMP_NUM_THREADS=4 ./build/tools/PATO/PATO ...
```

### Setting the number of simultaneous sequences

To reduce the memory footprint of PATO: high PerformAnce TriplexatOr, one can
set the maximum number of sequences that may processed simultaneously by the
triplex search algorithms. This is done by setting the `-cs` or `--chunk-size`
option to a value greater than 0 (128 by default). For instance, to process a
dataset in chunks of 32 sequences, execute:

```bash
$ ./build/tools/PATO/PATO --chunk-size 32 ...
```

To give an upper bound of the memory consumption of PATO: high PerformAnce
TriplexatOr, one can use the following formula:

$$
\begin{flalign}
& \text{mem}(cs, t, l) = cs \cdot \text{tpx}(l) + cs \cdot \text{len}(l) + \begin{cases}
t \cdot \text{len}(l), & \text{if } cs > t\\
cs \cdot \text{len}(l), & \text{if } cs \leq t
\end{cases} &
\end{flalign}
$$

where $cs$ is the number of simultaneous sequences, $t$ is the number of
threads, $l$ is the longest sequence in a given dataset, and $\text{tpx}(l)$ is
the size of the triplex features found in $l$.

It is possible to further reduce the memory usage of the application by
disabling the filtering of low-complex regions. This can be done by setting the
`-fr` or `--filter-repeats` option to `false`. In such a case, sequences should
be filtered before being passed to PATO: high PerformAnce TriplexatOr (the
[Ensembl genome browser](https://www.ensembl.org) provides filtered sequences)
and the formula becomes:

$$
\begin{flalign}
& \text{mem}(cs, t, l) = cs \cdot \text{tpx}(l) + cs \cdot \text{len}(l) &
\end{flalign}
$$

In general, one can't go wrong by setting the number of simultaneous sequences
to a value equal to the number of threads that PATO: high PerformAnce
TriplexatOr is going to use. However, if the sequences of a dataset are very
long, it may be necessary to reduce the number of simultaneous sequences to
avoid running out of memory.

If you are unsure about the number of simultaneous sequences to use, you can set
the `-cs` or `--chunk-size` option to 1. Although this may hurt parallelism by a
small amount, it will allow you to run PATO: high PerformAnce TriplexatOr on any
dataset without having to worry about the memory footprint of the application.

## Cite us

If you use PATO: high PerformAnce TriplexatOr in your research, please cite our
work using the following reference:

```bibtex
@article{amatria2023pato,
  title={PATO: genome-wide prediction of {lncRNA--DNA} triple helices},
  author={Amatria-Barral, I{\~n}aki and Gonz{\'a}lez-Dom{\'\i}nguez, Jorge and Touri{\~n}o, Juan},
  journal={Bioinformatics},
  volume={39},
  number={3},
  pages={btad134},
  year={2023}
}
```

## License

PATO: high PerformAnce TriplexatOr is free software and as such it is
distributed under the [MIT License](LICENSE).
