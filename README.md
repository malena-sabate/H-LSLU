# LSLU: Least Squares LU 

This is an inner product free Krylov Subspace Method.

This repository contains a MATLAB implementation of the **Least Squares LU (LSLU)** algorithm — a Krylov subspace method for solving large-scale rectangular linear inverse problems — along with a minimal toy example.

## Overview

**LSLU** is a modification of the Hessenberg iterative algorithm, designed for efficiently solving large-scale linear inverse problems. It is based on an LU factorization approach and operates in an **inner-product free** manner, making it particularly advantageous in high-performance computing environments and when using mixed precision arithmetic.

Moreover, **Hybrid LSLU**, incorporates **Tikhonov regularization** in a computationally efficient way.

### Key Features

- Solves large-scale rectangular linear inverse problems.
- Inner-product free: suitable for high-performance computing and mixed precision.
- Includes a hybrid variant with Tikhonov regularization.
- Competitive with existing iterative projection methods.

## Repository Contents

- `hybrid_LSLU.m` — Core implementation of the hybrid LSLU algorithm.
- `run_example.m` — A simple toy example demonstrating how to use the hybrid LSLU function.

> ⚠️ The toy example requires an external toolbox (see **Dependencies** below).

## Installation

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/LSLU.git
   cd LSLU

