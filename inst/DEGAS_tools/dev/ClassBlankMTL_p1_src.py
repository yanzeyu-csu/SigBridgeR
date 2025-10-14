"""
ClassBlankMTL Initialization Script (Part 1) - Source Version

This model is for: Single-cell data with classification labels, patient data without labels

Model Type: ClassBlank MTL
- Class: Single-cell data with classification labels
- Blank: Bulk data without labels (unlabeled)

This file will be processed by build_scripts.py
Utils functions from `mtl_utils.py` will be automatically inlined
"""

# INLINE_UTILS_HERE

# ***********************************************************************
# Set random seed
# ***********************************************************************

if sys.argv[10] != "NULL":
    np.random.seed(int(sys.argv[10]))
    tf.compat.v1.set_random_seed(int(sys.argv[10]))


# ***********************************************************************
# Data loading
# ***********************************************************************

data_folder = sys.argv[1]

# Load single-cell expression data and labels
Xsc = np.loadtxt(data_folder + "scExp.csv", delimiter=",", skiprows=1)
Ysc = np.loadtxt(data_folder + "scLab.csv", delimiter=",", skiprows=1)
Nsc = Ysc.shape[0]
Fsc = Xsc.shape[1]
Lsc = Ysc.shape[1]
idx_sc = np.arange(Nsc)
np.random.shuffle(idx_sc)

# Load patient expression data (no labels)
Xpat = np.loadtxt(data_folder + "patExp.csv", delimiter=",", skiprows=1)
Npat = Xpat.shape[0]
Fpat = Xpat.shape[1]
idx_pat = np.arange(Npat)
np.random.shuffle(idx_pat)


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
# Build network - Define placeholders
# ***********************************************************************

kprob = tf.placeholder(tf.float32)
xs = tf.placeholder(tf.float32, [None, Fsc])
ys_sc = tf.placeholder(tf.float32, [None, Lsc])
es = tf.placeholder(tf.float32, [None, Lsc])  # For autoencoder (if needed)
lsc = tf.placeholder(tf.int32, shape=())
lpat = tf.placeholder(tf.int32, shape=())

# ***********************************************************************

