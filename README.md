# Numerical validation for two related manuscripts

This repository contains R code used for the numerical analyses and figure generation associated with two related manuscripts on homeostatic mutation--selection balance and repair:

## 1. *Homeostatic Turnover, Repair, and the Demographic Form of Haldane's Principle*

**Ben Hardisty. *Homeostatic Turnover, Repair, and the Demographic Form of Haldane's Principle*.**

The code associated with this manuscript numerically validates the asymptotic and exact results derived for a homeostatically regulated, damage-structured population with reproduction-linked mutation and continuous repair.

## 2. *Repair-Mediated Path Filtering of Mutation Load in Multilocus Homeostatic Populations*

**Ben Hardisty. *Repair-Mediated Path Filtering of Mutation Load in Multilocus Homeostatic Populations*.**

The code associated with this manuscript numerically validates the multilocus path-filtering results, including the finite-mutation suppression ratio \(Q_k(U,\nu)\), convergence to the leading-order suppression factor \(E_k\), and the predicted independence of the leading suppression factor from the total number of loci for fixed target class \(k\).

The current canonical script for the path-filtering manuscript is:

`RepairMediatedPathFiltering_AllFigures_v4.R`

This script reproduces the main numerical figure and both supplementary numerical figures used in the manuscript.

Older scripts are retained for provenance but are superseded by the current canonical version.
---

## Requirements

The script requires R and the package:

```r
deSolve
