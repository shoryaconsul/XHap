import argparse

from CAECseq_haplo import train_caecseq

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filehead", help="Prefix of required files",
	type=str, required=True)
parser.add_argument("-p", "--ploidy", help="Ploidy of organism",
    default=2, type=int)
parser.add_argument("-n", "--num_expt", help="Number of datasets",
	default=1, type=int)
parser.add_argument("-a", "--algo_runs", help="Number of experimental runs per dataset",
	default=1, type=int)
parser.add_argument("--verbose", help="True for more verbose output",
    action='store_true', default=False)

args = parser.parse_args()
mec_expt = []
for i in range(args.num_expt):
    fhead = args.filehead + "_iter" + str(i+1)
    # fhead = args.filehead
    mec = train_caecseq(fhead + '/' + fhead, args.ploidy, args.algo_runs)

