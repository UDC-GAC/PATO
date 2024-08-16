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
