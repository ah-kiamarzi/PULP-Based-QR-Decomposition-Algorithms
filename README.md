### PULP-Based QR Decomposition Algorithms

This repository contains three implementations of the PULP-based QR decomposition algorithm: Given Rotation, Gram Schmidt, and Householder. The QR decomposition is a widely used linear algebra technique that can be used to solve a variety of problems, including least squares optimization, eigenvalue computation, and more.Each algorithm has been implemented in C and optimized for PULP-based architectures.

**### Background**

QR decomposition is a matrix factorization technique used to decompose a matrix into a product of an orthogonal matrix and an upper triangular matrix. It is commonly used in linear regression, signal processing, and data compression.
The three algorithms implemented in this repository differ in their approach to performing the QR decomposition. The Given Rotation algorithm uses a series of Givens rotations to transform the matrix into upper triangular form. The Gram Schmidt algorithm orthogonalizes the columns of the matrix using a modified Gram Schmidt process. The Householder algorithm uses a series of Householder reflections to transform the matrix into upper triangular form.

**### Dependencies**

The implementations in this repository require the PULP-SDK tool-chain to be installed. Please refer to the [PULP-RISCV-GNU-TOOLCHAIN](https://github.com/pulp-platform/pulp-riscv-gnu-toolchain) and [PULP_SDK](https://github.com/pulp-platform/pulp-sdk) for instructions on installing the PULP RISCV GNU TOOLCHAIN and the PULP SDK.

**### Building and Running**

To build each algorithm, navigate to the corresponding directory (Given_Rotation, Gram_Schmidt, or Householder) and use the following command:

```bash
make clean all run USE_CLUSTER=<UC> NUM_CORES=<#CORES>
```
This will generate an executable file. Where UC is a binary value, and #CORES is the number of cores that execute the code.

**### Performance**
The performance of each algorithm was evaluated using GVSoC platform. The results are shown below:

Algorithm
Matrix Size
Execution Time (ms)
Given Rotation


Gram Schmidt


Householder
















