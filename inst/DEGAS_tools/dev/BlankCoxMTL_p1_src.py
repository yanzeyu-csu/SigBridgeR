"""
# BlankCoxMTL Initialization Script (Part 1) - Source Version

This model is used for Cox survival analysis of patient data in the absence of labels in single-cell data.

Model Type: BlankCox MTL
- Blank: No single cell level labels
- Cox: The patient data is analyzed using the Cox proportional hazards model

This file will be processed by build_scripts.py
Utils functions from mtl_utils.py will be automatically inlined
"""

# INLINE_UTILS_HERE

# ***********************************************************************
# Seed
# ***********************************************************************

if sys.argv[10] != "NULL":
    np.random.seed(int(sys.argv[10]))
    tf.compat.v1.set_random_seed(int(sys.argv[10]))


# ***********************************************************************
# Data load
# ***********************************************************************

data_folder = sys.argv[1]

# Load single-cell expression data (unlabeled)
Xsc = np.loadtxt(data_folder + "scExp.csv", delimiter=",", skiprows=1)
Nsc = Xsc.shape[0]
Fsc = Xsc.shape[1]
idx_sc = np.arange(Nsc)
np.random.shuffle(idx_sc)

# Load patient expression data and survival data.
Xpat = np.loadtxt(data_folder + "patExp.csv", delimiter=",", skiprows=1)
Ypat = np.loadtxt(data_folder + "patLab.csv", delimiter=",", skiprows=1)
Npat = Xpat.shape[0]
Fpat = Xpat.shape[1]
Lpat = 1  # The output dimension of the Cox model is 1 (risk score).

idx_pat = np.arange(Npat)
np.random.shuffle(idx_pat)

survtime = Ypat[:, 0]
censor = Ypat[:, 1]


# ***********************************************************************
# Hyperparameters
# ***********************************************************************

train_steps = int(sys.argv[2])
scbatch_sz = int(sys.argv[3])
patbatch_sz = int(sys.argv[4])
hidden_feats = int(sys.argv[5])
do_prc = float(sys.argv[6])
lambda1 = float(sys.argv[7])
lambda2 = float(sys.argv[8])
lambda3 = float(sys.argv[9])


# ***********************************************************************
# Network placeholders
# ***********************************************************************

kprob = tf.placeholder(tf.float32)
xs = tf.placeholder(tf.float32, [None, Fsc])
r_pat = tf.placeholder(tf.float32, [None, None])  # Risk set matrix
c_pat = tf.placeholder(tf.float32, [None])        # Censoring indicator variable
ps = tf.placeholder(tf.float32, [None, Lpat])
lsc = tf.placeholder(tf.int32, shape=())
lpat = tf.placeholder(tf.int32, shape=())

# ***********************************************************************
