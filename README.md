# tsxCount
[![Build Status](https://travis-ci.org/mjoppich/tsxCount.svg?branch=master)](https://travis-ci.org/mjoppich/tsxCount)

## Introduction

tsxCount is an implementation of a HashMap specialised for k-mer counting.
K-mer counting is a task, where all k-long substrings of a string consisting of the letters A,T,C,G are counted.

Example for the string `ATCGAGTCAGTA` with all k=5-mers printed below:
```
ATCGAGTCAGTA
ATCGA
 TCGAG
  CGAGT
   GAGTC
    AGTCA
     GTCAG
      TCAGT
       CAGTA
```
The funny thing about k-mer counting is, that given the limited alphabet, no characters need to be stored, but each ATCG-character can be encoded in 2 bits - which reduces the memory requirement for storing these k-mers by 75%.

Unluckily, for genomics applications, the k-parameter wants to be chosen large (possibly k > 128). Unfortunately this means that 256bit are required to store the k-mer in, which does not fit any base integer class anymore, hence BigIntegers are required - inducing another burden, and requiring more steps in handling.

Furthermore, using a large array with one position for each k-mer requires too much allocated memory. The distribution of these k-mers, depending on the source of the strings, will be quite divisive: there will be a huge number of k-mers which occur only few times (once), but there will also be some k-mers which occur really frequently (hundred-thousands).
Thus, a more complex mechanism for encoding a) the k-mers and b) handling this distribution is required. [Marcais et al.](https://academic.oup.com/bioinformatics/article/27/6/764/234905) introduced such a data structure with their [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) implementation.
This structure has been re-implemented here, with some modifications.

## Scope

With tsxCount, multiple serialization techniques, which protect the increment of a specific k-mer are evaluated.
A specific implementation using hardware transactional memory via Intel TSX-NI (TSX) is compared with conventional approaches such as PTHREAD-mutexes or OpenMP-mutexes (OMP). In contrast to the PTHREAD mutexes, the OpenMP-mutexes allow lock hints, which here was set to [`omp_lock_hint_speculative`](https://www.openmp.org/spec-html/5.0/openmpsu155.html). Additionally, another lock-free approach is implemented making use of the Compare-And-Swap (CAS) atomic instruction, similar to jellyfish.

tsxCount is not meant to be a new k-mer counter. Instead, tsxCount is meant to compare serialization techniques a non-trivial application.

## Installation

tsxCount is cmake enabled. Hence, given that all dependencies are available (compiler, libgomp), the tsxCount build process should look similar to:

```
git clone https://github.com/mjoppich/tsxCount.git
cd tsxCount
mkdir build
cd build
CXX=clang++ cmake ..
make -j 4

./tsxCount --help
```


## Test Installation

A minimal working example is included in the data folder.

```
./tsxCount --input=../data/small_t7.1000.fastq --mode=OMP --check
```
The --check parameter ensures that the counts are correct (match with the `../data/small_t7.1000.fastq.14.count` file).

