# Introduction

Matrix-matrix multiplication is at the core of much of scientific computing. Despite many years of research, there are still open questions about the most efficient implementation mechanism. 

Lots of matrices that we encounter in the wild are sparse. That is, they contain many zero entries. In the limit of very large matrices, we do not have the memory to store them in a dense format. Even if we did, performing a matrix-matrix multiplication on this format would carry out a large amount of unnecessary computation. Moreover, many sparse matrices tend to exhibit some structure: large blocks of the matrix are completely empty (only zeroes), while other parts are quite populated. Many sparse matrices contain structured dense subblocks.
As an example, consider multiplying two n × n diagonal matrices. If we store them densely, this requires 2n^3 operations, whereas a sparse format would only require 2n operations.

This repository contains code written in C for optimised sparse matrix-matrix addition and multiplication for a 3rd year Assignment as part of my Masters Degree in Computer Science.

## Usage

Sample matrices are provided in the zip files.

```bash
cd GPU_Compute
...
```

More stuff here
