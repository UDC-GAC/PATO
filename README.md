# PATO: high PerformAnce TriplexatOr
[![AUR package](https://repology.org/badge/version-for-repo/aur/pato.svg)](https://repology.org/project/pato/versions)
[![compile and test gnu](https://github.com/amatria/pato/actions/workflows/compile-and-test-gnu.yml/badge.svg)](https://github.com/amatria/pato/actions/workflows/compile-and-test-gnu.yml)
[![compile and test clang](https://github.com/amatria/pato/actions/workflows/compile-and-test-clang.yml/badge.svg)](https://github.com/amatria/pato/actions/workflows/compile-and-test-clang.yml)

PATO: high PerformAnce TriplexatOr is a high performance tool for the fast and efficient detection of triple helices and triplex features in nucleotide sequences. PATO: high PerformAnce TriplexatOr, [Triplexator](https://github.com/Gurado/triplexator), and [TDF](https://github.com/CostaLab/reg-gen) are exact algorithms, and so PATO: high PerformAnce TriplexatOr functions nearly as a drop in replacement to accelerate the triplex analyses in multicore computers and achieves the same level of prediction accuracy as that of the state of the art.

## Version
Version v0.0.0.

## Requirements
To compile and execute PATO: high PerformAnce TriplexatOr, the following software is required:
* GNU Make
* C++ compiler (c++17 compliant and with support for OpenMP directives)

For instance, a valid combination of these tools may be: GNU Make v3.82, and GCC v9.3.0.

## Compilation
Download the source code from this repository, either use Git or download a copy from GitHub, and let GNU Make automatically compile PATO: high PerformAnce TriplexatOr for you:
```bash
user@host:/path/to/pato$ make gnu -j$(getconf _NPROCESSORS_ONLN)
```

## License
PATO: high PerformAnce TriplexatOr is free software and as such it is distributed under the [MIT License](LICENSE). However, PATO: high PerformAnce TriplexatOr makes use of several modules which are not original pieces of work. Therefore, their usage is subject to their corresponding [THIRDPARTLICENSE](THIRDPARTYLICENSES) and all rights are reserved to their authors.
