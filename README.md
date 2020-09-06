Algorithm
-------
info_incorporated GGM: Information-incorporated Gaussian Graphical Model

Maintainer
-------
Huangdi Yi (<huangdi.yi@yale.edu>)


Publication
-------
Yi H, Zhang Q, Lin C, Ma S (2020). Information-incorporated Gaussian Graphical Model for Gene Expression Data. Manuscript.


Usage
-------
1. ```info_incorporatedGGM.R```: main functions used in estimation of the Information-incorporated Gaussian Graphical Model

2. ```networks.R```: functions to generate Erdos-Renyi, scale-free, nearest-neighbor, and (positive and negative) banded network structures.

3. ```main.R```: an example under one simulation setting (p=100, n=300), containing
   * Generating the true precision matrix, the corresponding covariance matrix,
   * Generating the prior infomation
   * Estimation

Output
------
Four lists containing the estimation for each replicate under four different prior scenarios using the proposed info_incorporated GGM: ```Est.1```, ```Est.2```, ```Est.3```, and ```Est.4```.

Four lists containing the estimation for each replicate under four different prior scenarios using the infomation guided GGM (priorGGM, Jiang et al. (2016)): ```prior1```, ```prior2```, ```prior3```, and ```prior4```.
