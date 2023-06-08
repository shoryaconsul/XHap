import numpy as np
import logging
import typing
import shutil
from tqdm import tqdm  

import torch
from torch import nn
from torch import optim
from torch.utils.data import Dataset, DataLoader
from torch.nn import functional as F
from matplotlib import pyplot as plt


def save_ckp(state, checkpoint_path):
    """
    state: checkpoint we want to save
    checkpoint_path: path to save checkpoint
    """
    f_path = checkpoint_path  # Save path
    torch.save(state, f_path)


def load_ckp(checkpoint_path, model, optimizer):
    """
    checkpoint_path: path to save checkpoint
    model: model that we want to load checkpoint parameters into       
    optimizer: optimizer we defined in previous training
    """
    checkpoint = torch.load(checkpoint_path)
    # initialize state_dict from checkpoint to model
    model.load_state_dict(checkpoint['state_dict'])
    # initialize optimizer from checkpoint to optimizer
    optimizer.load_state_dict(checkpoint['optimizer'])

    # return model, optimizer, epoch value
    return model, optimizer, checkpoint['epoch']


class MyFilter(object):
    def __init__(self, level):
        self.__level = level

    def filter(self, logRecord):
        return logRecord.levelno <= self.__level


class SNVMatrixDataset(Dataset):
    def __init__(self, SNV_file, transform=None):
        """
        SNV_file: txt file containing SNV matrix
        """
        SNV_matrix_raw = np.loadtxt(SNV_file, dtype=int)
        self.SNV_matrix = SNV_matrix_raw[np.sum(SNV_matrix_raw != 0, axis=1) > 1]
	 
    def __len__(self):
        return np.shape(self.SNV_matrix)[0]
    
    def __getitem__(self, idx):
        SNV_row = torch.from_numpy(self.SNV_matrix[idx])
        SNV_row_onehot = F.one_hot(SNV_row, 5)[:,1:]
        SNV_row_onehot = SNV_row_onehot.type(torch.float32)
        SNV_row_onehot = SNV_row_onehot.transpose(1,0)
        return SNV_row_onehot[None,:], idx  # Shape is batch x 4 x numSNP


class ReadAE(nn.Module):
    def __init__(self, nSNP: int, latent_dim: int=None):
        super().__init__()
        self.nSNP = nSNP
        if latent_dim is None:
        	latent_dim = int(np.ceil(nSNP/4))  # Size of embedding

        self.encoder = nn.Sequential(
            nn.Conv2d(1, 32, (4,5), (4,1), (0, 2)),  # Padding changed from (1,2) to (0,2)
            nn.PReLU(),
            nn.Conv2d(32, 64, (1,5), (1,1), 'same'),
            nn.PReLU(),
            nn.Conv2d(64, 128, (1,3), (1,1), 'same'),
            nn.PReLU(),
            nn.Flatten(),
            )
        
        self.fc1 = nn.Linear(128*nSNP, latent_dim)
        self.fc2 = nn.Linear(latent_dim, 128*nSNP)
        self.act1 = nn.PReLU()

        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(128, 64, (1,3), (1,1), (0, 1)),
            nn.PReLU(),
            nn.ConvTranspose2d(64, 32, (1,5), (1,1), (0, 2)),
            nn.PReLU(),
            nn.ConvTranspose2d(32, 1, (4,5), (4,1), (0,2)),
            # nn.PReLU()
            )

    def forward(self, x):
        x_code = self.encoder(x)
        x_fc1 = self.fc1(x_code)
        x_flatten = self.act1(self.fc2(x_fc1))
        x_reshape = x_flatten.view(-1, 128, 1, self.nSNP)
        return x_fc1, self.decoder(x_reshape)


def learn_embed(dataset: Dataset, num_epoch: int,
				embed_dim: int=None,
				savefile: str = None,
				logger: logging.Logger = None,
				) -> ReadAE:
	"""
	dataset: 
		torch.utils.data.Dataset object
	num_epoch: 
		number of epochs to train network
	savefile:
		path to file for storing trained network weights
	logger:
		logging.Logger object to log progress
	"""
	# SNVdata = SNVMatrixDataset('Simulated_data/K3/cov15/sample4/simu_erro1_K3_cov15'\
	#                            '_l5000_iter_7_SNV_matrix.txt')

	# Setting up logging
	if logger is None:
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.INFO)

		# Create handlers
		c_handler = logging.StreamHandler()
		f_handler = logging.FileHandler('embed_train.log', mode='w')
		c_handler.setLevel(logging.WARNING)
		f_handler.setLevel(logging.INFO)
		f_handler.addFilter(MyFilter(logging.INFO))

		# Create formatters and add it to handlers
		c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
		f_format = logging.Formatter('%(asctime)s %(message)s')
		c_handler.setFormatter(c_format)
		f_handler.setFormatter(f_format)

		# Add handlers to the logger
		logger.addHandler(c_handler)
		logger.addHandler(f_handler)

	# Loading data
	SNVdata = dataset
	# num_epoch = 100
	nSNP = SNVdata[0][0].shape[2] # Number of SNVs
	num_read = len(SNVdata)  # Number of reads
	batch_size = int(np.ceil(num_read/20))

	dataloader = DataLoader(SNVdata, batch_size=batch_size,
	                        shuffle=True, num_workers=0)
	device = torch.device("cuda" if torch.cuda.is_available() else "cpu") #  use gpu if available

	embedAE = ReadAE(nSNP, embed_dim).to(device)  # Create and send model to device

	optimizer = optim.Adam(embedAE.parameters(), lr=1e-2)
	MSE = nn.MSELoss()
	train_loss_arr = []

	for epoch in tqdm(range(num_epoch)):
	    loss = 0
	    for batch_data, _ in dataloader:
	        optimizer.zero_grad()  # reset the gradients back to zero
	        _, recon = embedAE(batch_data) # compute reconstructions
	        train_loss = MSE(recon, batch_data)  # compute training reconstruction loss
	        train_loss.backward()  # compute accumulated gradients
	        optimizer.step()
	        loss += train_loss.item()  # add the mini-batch training loss to epoch loss
	    loss = loss / len(dataloader)  # compute the epoch training loss
	    train_loss_arr.append(loss)

	    # display the epoch training loss
	    logger.info("epoch : {}/{}, loss = {:.2f}".format(epoch + 1, num_epoch, loss))
	    # print("epoch : {}/{}, loss = {:.2f}".format(epoch + 1, num_epoch, loss))
	    if savefile and (epoch % 10 == 0):
	    	checkpoint = {
            'epoch': epoch + 1,
            'state_dict': embedAE.state_dict(),
            'optimizer': optimizer.state_dict(),
        	}
	    	save_ckp(checkpoint, savefile)
	return embedAE

if __name__ == "__main__":
	datapath = 'Simulated_data/diploid/cov15/sample1/simu_erro1_K2_cov5'\
				'_l5000_iter_1_SNV_matrix.txt'
	# datapath = 'Simulated_data/K3/cov15/sample4/simu_erro1_K3_cov15'\
	#                            '_l5000_iter_7_SNV_matrix.txt'
	SNVdata = SNVMatrixDataset(datapath)
	num_epoch = 100
	embedAE = learn_embed(SNVdata, num_epoch)
	print(embedAE)
