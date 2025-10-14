"""
ClassClassMTL Training Script (Part 3) - Source Version

This script defines the training loop and loss functions
For both single-cell and bulk classification tasks

This file will be processed by `build_scripts.py`
"""

# ***********************************************************************
# Prediction layer definitions
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

# Patient classification head
predict_pat = add_layer(
    layerF,
    hidden_feats,
    Lpat,
    activation_function=tf.nn.softmax,
    dropout_function=False,
    lambda1=lambda1,
)

# Autoencoder: SC to Patient mapping
layerae_sc2pat = add_layer(
    es,
    Lsc,
    hidden_feats,
    activation_function=tf.sigmoid,
    dropout_function=False,
    lambda1=lambda1,
)

predictae_sc2pat = add_layer(
    layerae_sc2pat,
    hidden_feats,
    Lpat,
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

# Patient classification loss
lossLabel2 = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(ys_pat - tf.slice(predict_pat, [lsc, 0], [lpat, Lpat])),
        reduction_indices=[1],
    )
)

# MMD loss for domain adaptation
lossMMD = mmd_loss(
    tf.slice(layerF, [0, 0], [lsc, hidden_feats]),
    tf.slice(layerF, [lsc, 0], [lpat, hidden_feats]),
)

# Constraint loss: SC samples should have uniform patient predictions
lossConstSCtoPT = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(tf.slice(predict_pat, [0, 0], [lsc, Lpat]) - (1.0 / Lpat)),
        reduction_indices=[1],
    )
)

# Constraint loss: Patient samples should have uniform SC predictions
lossConstPTtoSC = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(tf.slice(predict_sc, [lsc, 0], [lpat, Lsc]) - (1.0 / Lsc)),
        reduction_indices=[1],
    )
)

# Combined loss for main task
loss = (
    2 * lossLabel1
    + lambda2 * lossLabel2
    + lambda3 * lossMMD
    + lossConstSCtoPT
    + lossConstPTtoSC
)

# Autoencoder loss
lossae_sc2pat = tf.reduce_mean(
    tf.reduce_sum(tf.square(ps - predictae_sc2pat), reduction_indices=[1])
)

# Optimizers
train_step1 = tf.train.AdamOptimizer(learning_rate=0.01, epsilon=1e-3).minimize(loss)
train_step2 = tf.train.AdamOptimizer(learning_rate=0.01, epsilon=1e-3).minimize(
    lossae_sc2pat
)


# ***********************************************************************
# Training batch preparation function
# ***********************************************************************


def prepare_training_batch():
    """
    Prepare training batch with class balancing and validation

    Returns:
        tuple: (tensor_train, train_pat2) where
            - tensor_train: dict of training tensors
            - train_pat2: indices of selected patient samples
    """
    # Resample to balance classes
    train_sc = resample(50, Ysc, idx_sc)
    train_pat = resample(50, Ypat, idx_pat)

    # Shuffle
    np.random.shuffle(train_sc)
    np.random.shuffle(train_pat)

    # Select batch
    train_sc2 = train_sc[0:scbatch_sz]
    train_pat2 = train_pat[0:patbatch_sz]

    # Ensure each class has at least 2 samples (for stability)
    while (
        np.sum(np.sum(np.squeeze(Ypat[train_pat2, :]) > 0, axis=0) < 2) > 0
        or np.sum(np.sum(np.squeeze(Ysc[train_sc2, :]) > 0, axis=0) < 2) > 0
    ):
        train_sc = resample(50, Ysc, idx_sc)
        train_pat = resample(50, Ypat, idx_pat)
        np.random.shuffle(train_pat)
        np.random.shuffle(train_sc)
        train_sc2 = train_sc[0:scbatch_sz]
        train_pat2 = train_pat[0:patbatch_sz]

    # Apply Gamma resampling for data augmentation
    resampleGammaXYpat = resample_mixGamma(
        np.squeeze(Xpat[train_pat2, :]),
        np.squeeze(Ypat[train_pat2, :]),
        list(range(patbatch_sz)),
        patbatch_sz,
        Lpat,
    )

    resampleGammaXYsc = resample_mixGamma(
        np.squeeze(Xsc[train_sc2, :]),
        np.squeeze(Ysc[train_sc2, :]),
        list(range(scbatch_sz)),
        scbatch_sz,
        Lsc,
    )

    # Prepare training dictionary
    tensor_train = {
        xs: np.concatenate([resampleGammaXYsc[0], resampleGammaXYpat[0]]),
        ys_sc: resampleGammaXYsc[1],
        ys_pat: resampleGammaXYpat[1],
        lsc: resampleGammaXYsc[1].shape[0],
        lpat: resampleGammaXYpat[1].shape[0],
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
print(
    f"{'Step':<10} {'Total':<12} {'SC Loss':<12} {'Pat Loss':<12} {'MMD':<12} {'AE Loss':<12}"
)
print("-" * 70)

for i in range(train_steps + 1):
    # Main training step
    sess.run(train_step1, feed_dict=tensor_train)

    # Autoencoder training step
    sess.run(
        train_step2,
        feed_dict={
            es: sess.run(
                predict_sc,
                feed_dict={xs: np.squeeze(Xpat[train_pat2, :]), kprob: do_prc},
            ),
            ps: sess.run(
                predict_pat,
                feed_dict={xs: np.squeeze(Xpat[train_pat2, :]), kprob: do_prc},
            ),
            kprob: do_prc,
        },
    )

    # Log and resample every 50 steps
    if i % 50 == 0:
        # Calculate loss values
        loss_val = sess.run(loss, feed_dict=tensor_train)
        lossLabel1_val = sess.run(lossLabel1, feed_dict=tensor_train)
        lossLabel2_val = sess.run(lossLabel2, feed_dict=tensor_train)
        lossMMD_val = sess.run(lossMMD, feed_dict=tensor_train)

        # Calculate autoencoder loss
        lossae_val = sess.run(
            lossae_sc2pat,
            feed_dict={
                es: sess.run(
                    predict_sc,
                    feed_dict={xs: np.squeeze(Xpat[train_pat2, :]), kprob: do_prc},
                ),
                ps: sess.run(
                    predict_pat,
                    feed_dict={xs: np.squeeze(Xpat[train_pat2, :]), kprob: do_prc},
                ),
                kprob: do_prc,
            },
        )

        # Print progress
        print(
            f"{i:<10} {loss_val:<12.4f} {lossLabel1_val:<12.4f} {lossLabel2_val:<12.4f} {lossMMD_val:<12.4f} {lossae_val:<12.4f}"
        )

        # Prepare next batch (except for last iteration)
        if i < train_steps:
            tensor_train, train_pat2 = prepare_training_batch()


