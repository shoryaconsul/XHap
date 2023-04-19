from sklearn.datasets import make_blobs
from matplotlib import pyplot as plt

import torch
import time

class KernelKMeans():
    def __init__(self, n_clusters: int, max_iter: int = 100, tol: float = 1e-3,
                 verbose: bool = False, device: torch.cuda.device=torch.device("cpu")):
        self.n_clusters = n_clusters  # Number of clusters
        self.max_iter = max_iter  # Cap on number of iterations if convergence not achieved
        self.tol = tol
        self.verbose = verbose  # Set to True if debugging output is required
        self.device = device  # 'cpu' to run on CPU, 'cuda' to run on GPU 

    def compute_dist(self, kernel_matrix: torch.Tensor):
        """
        Function to compute distances between points and cluster centers from given
        kernel matrix for data points.

        Parameters
        ----------
        kernel_matrix: 2-D PyTorch Tensor of floats representing kernel matrix

        Returns
        -------
        2-D PyTorch Tensor of distances between data points and cluster centers

        """

        m = kernel_matrix.size(dim=0)
        dist_matrix = torch.zeros((m, self.n_clusters), dtype=kernel_matrix.dtype, device=self.device)
        self.within_dist = torch.zeros(self.n_clusters, device=self.device)

        M = torch.zeros_like(dist_matrix, device=self.device)  # One hot cluster encoding
        M[torch.arange(m), self.labels_] = 1
        Msum = M.T @ M
        zero_idx = (torch.diag(Msum) == 0)  # Aboiding division by zero
        if torch.any(zero_idx).item():
            Msum[zero_idx, zero_idx] = 1

        # Vectorization to compute distance between points and cluster means
        Msum_inv = torch.diag(1./torch.diag(Msum))
        K = kernel_matrix

        wdist = torch.diag(M.T @ K @ M)*torch.diag(Msum_inv)**2
        self.within_dist = wdist

        dist_matrix = torch.outer(torch.ones(m, dtype=kernel_matrix.dtype, device=self.device), wdist
                                ) - 2*K @ M @ Msum_inv
        dmin = torch.amin(dist_matrix[:, ~zero_idx])
        dmax = torch.amax(dist_matrix[:, ~zero_idx])
        dist_matrix[:, zero_idx] = dmin + (dmax - dmin
                                           )*torch.rand(m, torch.sum(zero_idx),
                                                        dtype=kernel_matrix.dtype,
                                                        device=self.device)

        return dist_matrix

    def fit(self, kernel_matrix: torch.Tensor):
        """
        This function fits kernel K-means using given kernel K-means.

        Parameters
        ----------
        kernel_matrix: 2-D PyTorch Tensor of floats representing kernel matrix

        """

        n_data = kernel_matrix.size(dim=0)  # Number of data points
        self.labels_ = torch.randint(self.n_clusters, size=(n_data,),
                                    device=self.device)  # Random initialization

        for it in range(self.max_iter):

            dist_matrix = self.compute_dist(kernel_matrix)
            labels_prev = torch.clone(self.labels_)
            self.labels_ = torch.argmin(dist_matrix, axis=1)

            n_labels_changed = torch.count_nonzero(self.labels_ - labels_prev)
            # print(it, n_labels_changed/n_data)
            if n_labels_changed/n_data < self.tol:
                if self.verbose:
                    print("Converged at iteration %d with $.3f labels"
                          "changing" %(it, n_labels_changed/n_data))
                break
        return self


if __name__ == '__main__':
    n_cl = 2
    
    X, y = make_blobs(n_samples=1000, centers=n_cl, random_state=0)
    
    X = torch.from_numpy(X)
    y = torch.from_numpy(y)

    km = KernelKMeans(n_clusters=n_cl, max_iter=1000, tol=1e-3, device='cuda')
    Xn = torch.diag(1./torch.linalg.norm(X, axis=1)) @ X
    kernel_matrix = Xn @ Xn.T

    kernel_matrix = kernel_matrix.to(torch.device('cuda')).float()
    print(kernel_matrix.dtype)
    t0 = time.time()
    km.fit(kernel_matrix)
    t1 = time.time()
    print('Time taken for clustering = %.3fs' %(t1-t0))
    
    labels = km.labels_.to('cpu').numpy()
    plt.figure()
    for i in range(n_cl):
        plt.scatter(X[labels==i, 0], X[labels==i, 1])

    plt.figure()
    for i in range(n_cl):
        plt.scatter(X[y==i, 0], X[y==i, 1])