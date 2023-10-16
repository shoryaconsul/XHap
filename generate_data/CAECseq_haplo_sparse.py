import numpy as np
from sklearn.cluster import KMeans
from Bio import SeqIO
from helper_caec import *
import os
import sys
import argparse
import time

import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras.layers import Layer, InputSpec


# This class is  adopoted from Guo X., Liu X., Zhu E., Yin J. 2017. Deep Clustering with Convolutional Autoencoders.
class ClusteringLayer(Layer):
    def __init__(self, n_clusters, weights = None, alpha = 1.0, **kwargs):
        if 'input_shape' not in kwargs and 'input_dim' in kwargs:
            kwargs['input_shape'] = (kwargs.pop('input_dim'),)
        super(ClusteringLayer, self).__init__(**kwargs)
        self.n_clusters = n_clusters
        self.alpha = alpha
        self.initial_weights = weights
        self.input_spec = InputSpec(ndim = 2)

    def build(self, input_shape):
        assert len(input_shape) == 2
        input_dim = int(input_shape[1])
        self.input_spec = InputSpec(dtype = K.floatx(), shape = (None, input_dim))
        self.clusters = self.add_weight(shape = (self.n_clusters, input_dim), initializer = 'glorot_uniform', name = 'clusters')
        if self.initial_weights is not None:
            self.set_weights(self.initial_weights)
            del self.initial_weights
        self.built = True

    def call(self, inputs, **kwargs):
        q = 1.0 / (1.0 + (K.sum(K.square(K.expand_dims(inputs, axis = 1) - self.clusters), axis = 2) / self.alpha))
        q **= (self.alpha + 1.0) / 2.0
        q = K.transpose(K.transpose(q) / K.sum(q, axis = 1))
        
        return q


def train_caecseq(zonename, ploidy, expt, sparse=False):
    def caecseq_chunk(SNVonehot, SNVmatrix):
        n_read, _, n_SNV, _ = SNVonehot.shape
        n_clusters = ploidy
        drop_prob = 0
        learning_rate = 10 ** -3
        lam = 0.1
        update_interval = 20 * 1
        n_pretrain_epoch = 100
        itermax = 2 * 10 ** 4
        batch_size = int(np.ceil(n_read / 20))
        len_code =  int(np.ceil(n_SNV / 4))
        kernel = [5, 5, 3]
        filters = [32, 64, 128]
        mec_record = []
        haplo_record = []
        for num_experiment in range(expt):
            print('Running experiment: {}'.format(num_experiment + 1))
            # convolutional auto-encoder
            Input = tf.keras.Input(shape = (4, n_SNV, 1),
                                name = 'Input')
            # convolutional layer - 1
            conv_1 = tf.keras.layers.Conv2D(filters = filters[0],
                                            kernel_size = [4, kernel[0]],
                                            strides = (4, 1),
                                            padding = 'same',
                                            name = 'conv_1')(Input)
            PRelu_1 = tf.keras.layers.PReLU(name = 'PRelu_1')(conv_1)
            drop_1 = tf.keras.layers.Dropout(drop_prob, name = 'drop_1')(PRelu_1)

            conv_2 = tf.keras.layers.Conv2D(filters = filters[1],
                                            kernel_size = [1, kernel[1]],
                                            strides = (1, 1),
                                            padding = 'same',
                                            name = 'conv_2')(drop_1)
            PRelu_2 = tf.keras.layers.PReLU(name = 'PRelu_2')(conv_2)
            drop_2 = tf.keras.layers.Dropout(drop_prob, name = 'drop_2')(PRelu_2)

            conv_3 = tf.keras.layers.Conv2D(filters = filters[2],
                                            kernel_size = [1, kernel[2]],
                                            strides = (1, 1),
                                            padding = 'same',
                                            name = 'conv_3')(drop_2)
            PRelu_3 = tf.keras.layers.PReLU(name = 'PRelu_3')(conv_3)
            drop_3 = tf.keras.layers.Dropout(drop_prob, name = 'drop_3')(PRelu_3)
            # flatten
            flatten = tf.keras.layers.Flatten(name = 'flatten')(drop_3)
            # code
            code = tf.keras.layers.Dense(units = len_code,
                                        name = 'code')(flatten)
            # dense
            dense_1 = tf.keras.layers.Dense(units = flatten.shape[1],
                                            name = 'dense_1')(code)
            PRelu_4 = tf.keras.layers.PReLU(name = 'PRelu_4')(dense_1)
            # reshape
            reshape = tf.keras.layers.Reshape((drop_3.shape[1], drop_3.shape[2], drop_3.shape[3]), 
                                            name = 'reshape')(PRelu_4)
            # transposed convolution layer
            convT_1 = tf.keras.layers.Conv2DTranspose(filters = filters[1],
                                                    kernel_size = [1, kernel[2]],
                                                    strides = (1, 1),
                                                    padding = 'same',
                                                    name = 'convT_1')(reshape)
            PRelu_5 = tf.keras.layers.PReLU(name = 'PRelu_5')(convT_1)
            drop_4 = tf.keras.layers.Dropout(drop_prob, name = 'drop_4')(PRelu_5)

            convT_2 = tf.keras.layers.Conv2DTranspose(filters = filters[0],
                                                    kernel_size = [1, kernel[1]],
                                                    strides = (1, 1),
                                                    padding = 'same',
                                                    name = 'convT_2')(drop_4)
            PRelu_6 = tf.keras.layers.PReLU(name = 'PRelu_6')(convT_2)
            drop_5 = tf.keras.layers.Dropout(drop_prob, name = 'drop_5')(PRelu_6)

            convT_3 = tf.keras.layers.Conv2DTranspose(filters = 1,
                                                    kernel_size = [4, kernel[0]],
                                                    strides = (4, 1),
                                                    padding = 'same',
                                                    name = 'convT_3')(drop_5)
            drop_6 = tf.keras.layers.Dropout(drop_prob, name = 'drop_6')(convT_3)

            # clustering layer
            clustering_layer = ClusteringLayer(n_clusters, name = 'clustering_layer')(code)

            model = tf.keras.Model(inputs = Input,
                                outputs = [drop_6, clustering_layer])
            model.compile(loss = ['mean_squared_error', 'kld'],
                        loss_weights = [1 - lam, lam],
                        optimizer = tf.keras.optimizers.Adam(lr = learning_rate))

            # pretrain
            pretrain = tf.keras.Model(Input, drop_6)
            pretrain.compile(loss = 'mean_squared_error',
                            optimizer = tf.keras.optimizers.Adam(lr = learning_rate))

            pretrain_history = pretrain.fit(x = SNVonehot,
                                            y = SNVonehot,
                                            batch_size = batch_size,
                                            epochs = n_pretrain_epoch,
                                            verbose =  None)

            # initialize cluster centers using k-means
            encoder = tf.keras.Model(Input, code)
            pretrain_mec = []
            pretrain_center = []
            for i in range(10):
                kmeans = KMeans(n_clusters = n_clusters, n_init = 30)
                y_pred = kmeans.fit_predict(encoder.predict(SNVonehot))
                k_means_haplotypes = origin2haplotype(y_pred, SNVmatrix, n_clusters)
                pretrain_mec.append(MEC(SNVmatrix, k_means_haplotypes))
                pretrain_center.append(kmeans.cluster_centers_)
            index = np.argmin(pretrain_mec)
            model.get_layer(name = 'clustering_layer').set_weights([pretrain_center[index]])

            index = 0
            mec = []

            for i in range(itermax):
                if i % update_interval == 0:
                    _, q = model.predict(SNVonehot)
                    p = target_distribution_haplo(q, SNVmatrix, n_clusters)
                    y_pred = q.argmax(1)
                    haplotypes = origin2haplotype(y_pred, SNVmatrix, n_clusters)
                    mec.append(MEC(SNVmatrix, haplotypes))
                    if len(mec) > 1 and mec[-1] == mec[-2]:
                        break
                    
                if (index + 1) * batch_size > n_read:
                    if index * batch_size == n_read:
                        loss = model.train_on_batch(x = SNVonehot[index * batch_size - 1::],
                                                    y = [SNVonehot[index * batch_size - 1::], p[index * batch_size - 1::]])
                    else:
                        loss = model.train_on_batch(x = SNVonehot[index * batch_size::],
                                                    y = [SNVonehot[index * batch_size::], p[index * batch_size::]])
                    index = 0
                else:
                    loss = model.train_on_batch(x = SNVonehot[index * batch_size:(index + 1) * batch_size],
                                                y = [SNVonehot[index * batch_size:(index + 1) * batch_size],
                                                    p[index * batch_size:(index + 1) * batch_size]])
                    index += 1

            # correction
            pre_mec = 0
            mec = MEC(SNVmatrix, haplotypes)
            count = 0
            while mec != pre_mec:
                index = []
                for i in range(SNVmatrix.shape[0]):
                    dis = np.zeros((haplotypes.shape[0]))
                    for j in range(haplotypes.shape[0]):
                        dis[j] = hamming_distance(SNVmatrix[i, :], haplotypes[j, :])
                    index.append(np.argmin(dis))

                new_haplo = np.zeros((haplotypes.shape))
                for i in range(haplotypes.shape[0]):
                    reads_single = SNVmatrix[np.array(index) == i, :]
                    single_sta = np.zeros((n_SNV, 4))
                    if len(reads_single)!= 0:
                        single_sta = ACGT_count(reads_single)
                    new_haplo[i, :] = np.argmax(single_sta, axis = 1) + 1
                    uncov_pos = np.where(np.sum(single_sta, axis=1) == 0)[0]
                    if len(uncov_pos) != 0:
                        new_haplo[i, uncov_pos] = 0   

                haplotypes = new_haplo.copy()
                pre_mec = mec
                mec = MEC(SNVmatrix, haplotypes)
                count += 1
                
            mec_record.append(mec)
            haplo_record.append(new_haplo)
            
        haplotypes = haplo_record[np.argmin(mec_record)]
        return haplotypes

    zone_name = zonename # zone name
    SNVmatrix_name = zone_name + '_SNV_matrix.txt' # SNP matrix file name
    if not sparse:
        SNVmatrix, SNVonehot = import_SNV(SNVmatrix_name)
        haplotypes = caecseq_chunk(SNVonehot, SNVmatrix)
        print('Best MEC score: %.3f' %(MEC(SNVmatrix, haplotypes)))
    else:
        SNVmatrix_full, SNVonehot_full = import_sparseSNV(SNVmatrix_name)

        chunk_size = 250
        chunk_overlap = int(0.2*250)
        chunk_starts = list(range(0, SNVmatrix_full.shape[1],
                                  chunk_size - chunk_overlap))
        
        haplotypes = np.zeros((ploidy, SNVmatrix_full.shape[1]))  # Complete hapltoypes
        print("Number of chunks: %d" %len(chunk_starts))

        recon_end = 0
        t_0 = time.time()
        for i, cs in enumerate(chunk_starts):
            print("Running CAECSeq on chunk starting at SNP %d." %cs)
            SNVmatrix_chunk = SNVmatrix_full[:, cs:cs+chunk_size]
            chunk_rows = SNVmatrix_chunk.sum(axis=1).nonzero()[0]
            chunk_haplotypes = caecseq_chunk(
                SNVonehot_full[chunk_rows, :, cs:cs+chunk_size, :],
                SNVmatrix_chunk[chunk_rows]
                )
            if i == 0:
                haplotypes[:, cs:cs+chunk_size] = chunk_haplotypes
            else:
                prev_hap = haplotypes[:, cs:recon_end]
                match = best_match(prev_hap, chunk_haplotypes[:, :recon_end - cs])
                haplotypes[:, cs:cs + chunk_haplotypes.shape[1]] = chunk_haplotypes[match, :]
            recon_end = cs + chunk_size
            print("Finished running on %d chunks in %d seconds." %(i+1, time.time() - t_0))

    # Print reconstructed haplotypes to file
    with open(zone_name + '_Reconstructed_Strains.txt', 'w') as f:
        for i in range(haplotypes.shape[0]):
            f.write('Haplotype ' + str(i + 1) + '\n')
            for j in range(haplotypes.shape[1]):
                if haplotypes[i, j] == 1:
                    f.write('A')
                elif haplotypes[i, j] == 2:
                    f.write('C')
                elif haplotypes[i, j] == 3:
                    f.write('G')
                elif haplotypes[i, j] == 4:
                    f.write('T')
                elif haplotypes[i, j] == 0:
                    f.write('-')
            f.write('\n')

    print("CAECSeq finished!")
    # print("----------------------------------")
    # return np.amin(mec_record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("zonename", help="Prefix of SNV matrix file")
    parser.add_argument("-k", "--ploidy", help="Ploidy of haplotypes",
            type=int, required=True)
    parser.add_argument("-n", "--expt", help="Number of runs",
            type=int, default=1)
    parser.add_argument("-s", "--sparse", help="Use sparse reprsentation",
                        action='store_true', default=False)
    parser.add_argument("-p", "--gpu", help="GPU to use",
            default="-1")
    args = parser.parse_args()
    # print(args)

    os.environ["CUDA_DEVICE_ORDER"] = 'PCI_BUS_ID'
    os.environ['CUDA_VISIBLE_DEVICES']= args.gpu

    start_time = time.time()
    train_caecseq(args.zonename, args.ploidy, args.expt, args.sparse)
    print('Total time taken: %d seconds' %(time.time() - start_time))