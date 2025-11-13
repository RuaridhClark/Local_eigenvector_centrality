# Local Eigenvector Centrality

This repository contains code and example datasets for the publication, **"A Local Eigenvector Centrality"**.

The paper introduces a **new centrality metric** that extends traditional eigenvector centrality by incorporating both **global and local influence**.

---

## Repository Structure

| Folder / File                     | Description |
|----------------------------------|-------------|
| `High School contacts/`           | Example contact networks for a high school (MATLAB format) |
| `Primary School contacts/`        | Example contact networks for a primary school (MATLAB format) |
| `Road network/`                   | Example road network dataset (Google Colab notebook) |
| `local_eigenvector_centrality.m`  | MATLAB function to compute local eigenvector centrality |

**High School contact networks** from:  
  J. Stehlé, N. Voirin, A. Barrat, C. Cattuto, L. Isella, J.-F. Pinton, W. Q. Marco, V. den Broeck,  
  C. Régis, P. L. Bruno, and Vanhems. *High-resolution measurements of face-to-face contact patterns in a primary school*.  
  **PLOS ONE, 6:e23176, 2011.**  

**Primary School contact networks** from:  
  R. Mastrandrea, J. Fournet, and A. Barrat. *Contact patterns in a high school: A comparison between data collected using wearable sensors, contact diaries and friendship surveys.*  
  **PLOS ONE, 10:1–26, 2015.**  

The road network example call Open Street Maps and is implemented in a Jupyter notebook (`.ipynb`) for Google Colab

---

## MATLAB Function: `local_eigenvector_centrality.m`

Computes the **local eigenvector centrality** of a network adjacency matrix.

### Usage

```matlab
% Basic usage
centrality = local_eigenvector_centrality(A);

% Advanced usage with plotting and options
[centrality, details] = local_eigenvector_centrality(A, X, plotall, Imax);
```

### Inputs
```
A — Adjacency matrix (sparse or full)

X — (optional) Node coordinates [n x 2] for plotting

plotall — (optional) If true, generate plots visualising component eigenvectors and comparison with Bonacich's eigenvector centrality

Imax — (optional) Identified prominent eigengap
```
### Outputs
```
centrality — Local eigenvector centrality values

details    - Struct containing eigenvalues, eigenvectors, eigengap, and global eigenvector centrality for analysis/plotting
```
