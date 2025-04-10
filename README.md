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

1. To clone this repository:

   ```bash
   git clone https://github.com/yourusername/LSLU.git
   cd LSLU
   
### Software language

       MATLAB 9.14 (R2023a)
       For those without access to MATLAB, Octave provides an alternative platform.  
       Note that these codes have not been tested in Octave. 
       
### Dependencies

To run the toy example, you will need the following MATLAB toolboxes:

    "IR tools: A MATLAB Package of Iterative Regularization"
    by Silvia Gazzola, Per Christian Hansen and James G. Nagy
    https://github.com/jnagy1/IRtools.git

Make sure to download this toolbox, place the files in your MATLAB path, and call addpath as needed.

## Contributors
        Ariana N. Brown, 
        Department of Mathematics, Emory University

        Julianne Chung, 
        Department of Mathematics, Emory University
        
        James G. Nagy, 
        Department of Mathematics, Emory University
      
        Malena Sabaté Landman, 
        Mathematical Institute, University of Oxford
	
## Licensing

If you use this codes, you *must* cite the original authors:

       [1] Brown et al. "Inner Product Free Krylov Methods for Large-Scale Inverse Problems". 2025.

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

