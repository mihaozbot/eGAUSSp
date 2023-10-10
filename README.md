# Cluster-Volume-Based Merging Approach for Incrementally Evolving Fuzzy Gaussian Clustering—eGAUSS+

### Citation

Please cite the original paper as:
```
I. Škrjanc, “Cluster-Volume-Based Merging Approach for Incrementally Evolving Fuzzy Gaussian Clustering-eGAUSS+,” IEEE Transactions on Fuzzy Systems, vol. 28, no. 9, pp. 2222–2231, Sep. 2020, doi: 10.1109/TFUZZ.2019.2931874.
```
### Abstract

In this article, a new dynamic merging approach for incrementally evolving clustering is presented. The clusters are learned incrementally online from streams of data. The merging criterion is based on the comparison between the sum of volumes of two clusters and the expected volume of the newly merged cluster. The new merged cluster is formed using weighted averaging of cluster centers and the calculation of the joint covariance matrix. The proposed evolving algorithm eGAUSS+ is easy to implement, works on high-dimensional data sets, and produces reliable clusters.

### Repository Contents

- `Z2D.mat`: Contains pre-generated data
- `display_clusters.m`: MATLAB function to display clusters.
- `eGAUSSp.m`: Main script to run the eGAUSS+ algorithm.
- `eGAUSSp.mat`: Contains computed parameters

### Implementation

The code in this repository is the MATLAB 2021b implementation of the algorithm presented in the paper.

### Results
[egaussp.pdf](https://github.com/mihaozbot/eGAUSSp/files/12844826/egaussp.pdf)
