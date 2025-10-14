"""
BlankCoxMTL Training Script (Part 3) - Source Version

This script defines the training loop and loss functions
For single-cell without labels and patient survival analysis

This file will be processed by build_scripts.py
"""

# ***********************************************************************
# Prediction layer definition
# ***********************************************************************

predict_pat = add_layer(
    layerF,
    hidden_feats,
    Lpat,
    activation_function=tf.sigmoid,
    dropout_function=False,
    lambda1=lambda1,
)


# ***********************************************************************
# Loss functions
# ***********************************************************************

# Cox proportional hazards loss (negative log partial likelihood)
lossLabel2 = -tf.reduce_mean(
    (
        tf.squeeze(tf.slice(predict_pat, [lsc, 0], [lpat, Lpat]))
        - tf.log(
            tf.reduce_sum(
                tf.exp(tf.squeeze(tf.slice(predict_pat, [lsc, 0], [lpat, Lpat])))
                * r_pat,
                1,
            )
        )
    )
    * c_pat
)

# Maximum Mean Discrepancy loss for domain adaptation
lossMMD = mmd_loss(
    tf.slice(layerF, [0, 0], [lsc, hidden_feats]),
    tf.slice(layerF, [lsc, 0], [lpat, hidden_feats]),
)

# Constraint loss: SC predictions should be uniform
lossConstSCtoPT = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(tf.slice(predict_pat, [0, 0], [lsc, Lpat]) - (1.0 / 2.0)),
        reduction_indices=[1],
    )
)

# Combined loss
loss = 2 * lossLabel2 + lambda3 * lossMMD + lossConstSCtoPT

# Optimizer
train_step1 = tf.train.AdamOptimizer(learning_rate=0.01, epsilon=1e-3).minimize(loss)


# ***********************************************************************
# Training batch preparation function
# ***********************************************************************


def prepare_training_batch():
    """
    Prepare training batch with proper data shuffling

    Returns:
        dict: Training tensor dictionary containing:
            - xs: Concatenated SC and patient features
            - r_pat: Risk set matrix for Cox model
            - c_pat: Censoring indicators
            - lsc: Number of SC samples
            - lpat: Number of patient samples
            - kprob: Dropout probability
    """
    # Shuffle indices
    np.random.shuffle(idx_sc)
    np.random.shuffle(idx_pat)

    # Select batch samples
    train_sc2 = idx_sc[0:scbatch_sz]
    train_pat2 = idx_pat[0:patbatch_sz]

    # Prepare training dictionary
    tensor_train = {
        xs: np.concatenate(
            [np.squeeze(Xsc[train_sc2,]), np.squeeze(Xpat[train_pat2, :])]
        ),
        r_pat: Rmatrix(survtime[train_pat2]),
        c_pat: np.squeeze(censor[train_pat2]),
        lsc: len(train_sc2),
        lpat: len(train_pat2),
        kprob: do_prc,
    }

    return tensor_train


# ***********************************************************************
# Training main loop
# ***********************************************************************

# Initialize first batch
tensor_train = prepare_training_batch()

# Initialize TensorFlow session
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

print("Starting training...")
print(f"{'Step':<10} {'Total Loss':<15} {'Cox Loss':<15} {'MMD Loss':<15}")
print("-" * 55)

for i in range(train_steps + 1):
    # Training step
    sess.run(train_step1, feed_dict=tensor_train)

    # Log and resample every 50 steps
    if i % 50 == 0:
        # Calculate loss values
        loss_val = sess.run(loss, feed_dict=tensor_train)
        lossLabel2_val = sess.run(lossLabel2, feed_dict=tensor_train)
        lossMMD_val = sess.run(lossMMD, feed_dict=tensor_train)

        # Print progress
        print(f"{i:<10} {loss_val:<15.4f} {lossLabel2_val:<15.4f} {lossMMD_val:<15.4f}")

        # Prepare next batch (except for last iteration)
        if i < train_steps:
            tensor_train = prepare_training_batch()


