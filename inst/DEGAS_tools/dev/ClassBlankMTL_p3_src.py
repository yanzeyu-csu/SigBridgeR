"""
ClassBlankMTL Training Script (Part 3) - Source Version

This script defines the training loop and loss functions
For single-cell classification and unlabeled bulk data

This file will be processed by build_scripts.py
"""

# ***********************************************************************
# Prediction layer definition
# ***********************************************************************

# Single-cell classification head
predict_sc = add_layer(
    layerF,
    hidden_feats,
    Lsc,
    activation_function=tf.nn.softmax,
    dropout_function=False,
    lambda1=lambda1,
)


# ***********************************************************************
# Loss functions
# ***********************************************************************

# Single-cell classification loss
lossLabel1 = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(ys_sc - tf.slice(predict_sc, [0, 0], [lsc, Lsc])),
        reduction_indices=[1],
    )
)

# MMD loss for domain adaptation
lossMMD = mmd_loss(
    tf.slice(layerF, [0, 0], [lsc, hidden_feats]),
    tf.slice(layerF, [lsc, 0], [lpat, hidden_feats]),
)

# Constraint loss: Patient samples should have uniform SC predictions
lossConstPTtoSC = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(tf.slice(predict_sc, [lsc, 0], [lpat, Lsc]) - (1.0 / Lsc)),
        reduction_indices=[1],
    )
)

# Combined loss
loss = 2 * lossLabel1 + lambda3 * lossMMD + lossConstPTtoSC

# Optimizer
train_step1 = tf.train.AdamOptimizer(learning_rate=0.01, epsilon=1e-3).minimize(loss)


# ***********************************************************************
# Training batch preparation function
# ***********************************************************************

def prepare_training_batch():
    """
    Prepare training batch with SC class balancing
    
    Returns:
        tuple: (tensor_train, train_pat2) where
            - tensor_train: dict of training tensors
            - train_pat2: indices of selected patient samples (for compatibility)
    """
    # Resample SC data to balance classes
    train_sc = resample(50, Ysc, idx_sc)
    np.random.shuffle(train_sc)
    
    # Select batch
    train_sc2 = train_sc[0:scbatch_sz]
    train_pat2 = idx_pat[0:patbatch_sz]
    
    # Ensure each SC class has at least 2 samples
    while np.sum(np.sum(np.squeeze(Ysc[train_sc2, :]) > 0, axis=0) < 2) > 0:
        train_sc = resample(50, Ysc, idx_sc)
        np.random.shuffle(train_sc)
        train_sc2 = train_sc[0:scbatch_sz]
        
    
    # Apply Gamma resampling for SC data augmentation
    resampleGammaXYsc = resample_mixGamma(
        np.squeeze(Xsc[train_sc2, :]),
        np.squeeze(Ysc[train_sc2, :]),
        list(range(scbatch_sz)),
        scbatch_sz,
        Lsc,
    )
    
    # Prepare training dictionary
    tensor_train = {
        xs: np.concatenate([resampleGammaXYsc[0], np.squeeze(Xpat[train_pat2,])]),
        ys_sc: resampleGammaXYsc[1],
        lsc: resampleGammaXYsc[1].shape[0],
        lpat: len(train_pat2),
        kprob: do_prc,
    }
    
    return tensor_train, train_pat2


# ***********************************************************************
# Training main loop
# ***********************************************************************

# Initialize first batch
tensor_train, train_pat2 = prepare_training_batch()

# Initialize TensorFlow session
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

print("Starting training...")
print(f"{'Step':<10} {'Total Loss':<15} {'SC Loss':<15} {'MMD Loss':<15}")
print("-" * 55)

for i in range(train_steps + 1):
    # Training step
    sess.run(train_step1, feed_dict=tensor_train)
    
    # Log and resample every 50 steps
    if i % 50 == 0:
        # Calculate loss values
        loss_val = sess.run(loss, feed_dict=tensor_train)
        lossLabel1_val = sess.run(lossLabel1, feed_dict=tensor_train)
        lossMMD_val = sess.run(lossMMD, feed_dict=tensor_train)
        
        # Print progress
        print(f"{i:<10} {loss_val:<15.6f} {lossLabel1_val:<15.6f} {lossMMD_val:<15.6f}")
        
        # Prepare next batch (except for last iteration)
        if i < train_steps:
            tensor_train, train_pat2 = prepare_training_batch()


