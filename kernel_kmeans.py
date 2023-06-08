import numpy as np
from numpy import random as rn
from nptyping import NDArray
from typing import Any
import time


class KernelKMeans():
    def __init__(self, n_clusters: int, max_iter: int = 100, tol: float = 1e-3,
                 verbose: bool = False):
        self.n_clusters = n_clusters  # Number of clusters
        self.max_iter = max_iter  # Cap on number of iterations if convergence not achieved
        self.tol = tol
        self.verbose = verbose  # Set to True if debugging output is required

    def compute_dist(self, kernel_matrix: NDArray[(Any, Any), float]):
        m = kernel_matrix.shape[0]
        dist_matrix = np.zeros((m, self.n_clusters))
        self.within_dist = np.zeros(self.n_clusters)

        M = np.zeros_like(dist_matrix)  # One hot cluster encoding
        M[np.arange(m), self.labels_] = 1
        Msum = M.T @ M
        zero_idx = (np.diag(Msum) == 0)  # Aboiding division by zero
        if np.any(zero_idx):
            Msum[zero_idx, zero_idx] = 1

        # Vectorization to compute distance between points and cluster means
        Msum_inv = np.diag(1./np.diag(Msum))
        K = kernel_matrix
        wdist = np.diag(M.T @ K @ M)*np.diag(Msum_inv)**2
        self.within_dist = wdist

        dist_matrix = np.outer(np.ones(m), wdist) - 2*K @ M @ Msum_inv
        dmin = np.amin(dist_matrix[:, ~zero_idx])
        dmax = np.amax(dist_matrix[:, ~zero_idx])
        dist_matrix[:, zero_idx] = dmin + (dmax - dmin
                                           )*rn.rand(m, np.sum(zero_idx))

        return dist_matrix

    def fit(self, kernel_matrix: NDArray[(Any, Any), float]):
        n_data = kernel_matrix.shape[0]  # Number of data points
        self.labels_ = rn.randint(self.n_clusters, size=n_data)  # Random initialization

        for it in range(self.max_iter):

            dist_matrix = self.compute_dist(kernel_matrix)
            labels_prev = np.copy(self.labels_)
            self.labels_ = np.argmin(dist_matrix, axis=1)

            n_labels_changed = np.count_nonzero(self.labels_ - labels_prev)
            if n_labels_changed/n_data < self.tol:
                if self.verbose:
                    print("Converged at iteration %d with $.3f labels"
                          "changing" %(it, n_labels_changed/n_data))
                break
        return self


if __name__ == '__main__':
    from sklearn.datasets import make_blobs
    from matplotlib import pyplot as plt
    n_cl = 2
    
    X, y = make_blobs(n_samples=1000, centers=n_cl, random_state=0)
    km = KernelKMeans(n_clusters=5, max_iter=1000, tol=1e-3)
    Xn = np.diag(1./np.linalg.norm(X, axis=1)) @ X
    kernel_matrix = Xn @ Xn.T
    # kernel_matrix = np.zeros((X.shape[0], X.shape[0]))
    # for i in range(X.shape[0]):
    #     for j in range(X.shape[0]):
    #         kernel_matrix[i, j] = np.exp(-1*np.linalg.norm(X[i]-X[j])**2)
    
    t0 = time.time()
    km.fit(kernel_matrix)
    t1 = time.time()
    print('Time taken for clustering = %.3fs' %(t1-t0))
    
    plt.figure()
    for i in range(n_cl):
        plt.scatter(X[km.labels_==i, 0], X[km.labels_==i, 1])

    plt.figure()
    for i in range(n_cl):
        plt.scatter(X[y==i, 0], X[y==i, 1])