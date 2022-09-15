import numpy as np
from numpy import random as rn
from nptyping import NDArray
from typing import Any


class KernelKMeans():
    def __init__(self, n_clusters: int, max_iter: int = 100, tol: float = 1e-3,
                 verbose: bool = False):
        self.n_clusters = n_clusters  # Number of clusters
        self.max_iter = max_iter  # Cap on number of iterations if convergence not achieved
        self.tol = tol
        self.verbose = verbose  # Set to True if debugging output is required

    def compute_dist(self, kernel_matrix: NDArray[(Any, Any), float]):
        dist_matrix = np.zeros((kernel_matrix.shape[0], self.n_clusters))
        self.within_dist = np.zeros(self.n_clusters)

        for j in range(self.n_clusters):
            mask = (self.labels_ == j)
            # print(j, np.sum(self.labels_==j))
            if np.sum(mask) == 0:
                # raise ValueError("Empty cluster. Reinitializing cluster labels")
                dist_matrix[:, j] = ((np.argmax(dist_matrix)
                                     - np.argmin(dist_matrix))*rn.rand(kernel_matrix.shape[0])
                                     + np.argmin(dist_matrix))
            else:
                K_mask = kernel_matrix[mask][:, mask]
                dist_j = np.sum(K_mask) / np.sum(mask)**2
                self.within_dist[j] = dist_j
                dist_matrix[:, j] = (dist_matrix[:, j] + dist_j
                                     - 2*np.sum(kernel_matrix[:, mask], axis=1)/np.sum(mask))
        return dist_matrix

    def fit(self, kernel_matrix: NDArray[(Any, Any), float]):
        n_data = kernel_matrix.shape[0]  # Number of data points
        self.labels_ = rn.randint(self.n_clusters, size=n_data)  # Random initialization

        for it in range(self.max_iter):
            dist_matrix = self.compute_dist(kernel_matrix)
            labels_prev = np.copy(self.labels_)
            self.labels_ = np.argmin(dist_matrix, axis=1)

            n_labels_changed = np.count_nonzero(self.labels_ - labels_prev)
            # print(it, n_labels_changed/n_data)
            if n_labels_changed/n_data < self.tol:
                if self.verbose:
                    print("Converged at iteration %d with $.3f labels"
                          "changing" %(it, n_labels_changed/n_data))
                break
        return self


if __name__ == '__main__':
    from sklearn.datasets import make_blobs
    from matplotlib import pyplot as plt

    X, y = make_blobs(n_samples=1000, centers=5, random_state=0)
    km = KernelKMeans(n_clusters=5, max_iter=100, tol=1e-3)
    kernel_matrix = X @ X.T
    # kernel_matrix = np.zeros((X.shape[0], X.shape[0]))
    # for i in range(X.shape[0]):
    #     for j in range(X.shape[0]):
    #         kernel_matrix[i, j] = np.exp(-1*np.linalg.norm(X[i]-X[j])**2)
    
    km.fit(kernel_matrix)
    
    plt.figure()
    for i in range(5):
        plt.scatter(X[km.labels_==i, 0], X[km.labels_==i, 1])

    plt.figure()
    for i in range(5):
        plt.scatter(X[y==i, 0], X[y==i, 1])