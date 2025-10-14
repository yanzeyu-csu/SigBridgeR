"""
# BlankClassMTL Training Script (Part 3) - Source Version

This script defines the training loop and loss function, applicable to classification tasks with unlabeled single-cell data and labeled bulk data.

This file will be processed by `build_scripts.py`
"""

# ***********************************************************************
# Prediction Layer
# ***********************************************************************

predict_pat = add_layer(
    layerF,
    hidden_feats,
    Lpat,
    activation_function=tf.nn.softmax,
    dropout_function=False,
    lambda1=lambda1,
)


# ***********************************************************************
# Loss function
# ***********************************************************************

lossLabel2 = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(ys_pat - tf.slice(predict_pat, [lsc, 0], [lpat, Lpat])),
        reduction_indices=[1],
    )
)

lossMMD = mmd_loss(
    tf.slice(layerF, [0, 0], [lsc, hidden_feats]),
    tf.slice(layerF, [lsc, 0], [lpat, hidden_feats]),
)

lossConstSCtoPT = tf.reduce_mean(
    tf.reduce_sum(
        tf.square(tf.slice(predict_pat, [0, 0], [lsc, Lpat]) - (1.0 / Lpat)),
        reduction_indices=[1],
    )
)

loss = 2 * lossLabel2 + lambda3 * lossMMD + lossConstSCtoPT

train_step1 = tf.train.AdamOptimizer(learning_rate=0.01, epsilon=1e-3).minimize(loss)


# ***********************************************************************
# Training batch preparation function.
# ***********************************************************************

def prepare_training_batch():
    """
    Prepare training batches, including data resampling and validation.
    """
    train_pat = resample(50, Ypat, idx_pat)
    np.random.shuffle(train_pat)
    
    # Ensure that each class has at least two samples.
    while np.sum(np.sum(Ypat[train_pat[0:patbatch_sz], :], axis=0) < 2) > 0:
        np.random.shuffle(train_pat)
    
    np.random.shuffle(idx_sc)
    train_sc2 = idx_sc[0:scbatch_sz]
    train_pat2 = train_pat[0:patbatch_sz]
    
    resampleGammaXYpat = resample_mixGamma(
        np.squeeze(Xpat[train_pat2, :]),
        np.squeeze(Ypat[train_pat2, :]),
        list(range(patbatch_sz)),
        patbatch_sz,
        Lpat,
    )
    
    tensor_train = {
        xs: np.concatenate([np.squeeze(Xsc[train_sc2,]), resampleGammaXYpat[0]]),
        ys_pat: resampleGammaXYpat[1],
        lsc: len(train_sc2),
        lpat: resampleGammaXYpat[1].shape[0],
        kprob: do_prc,
    }
    
    return tensor_train


# ***********************************************************************
# Main training loop.
# ***********************************************************************

tensor_train = prepare_training_batch()
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

print("Starting training...")
print(f"{'Step':<10} {'Total Loss':<15} {'Label Loss':<15} {'MMD Loss':<15}")
print("-" * 55)

for i in range(train_steps + 1):
    sess.run(train_step1, feed_dict=tensor_train)
    
    if i % 50 == 0:
        loss_val = sess.run(loss, feed_dict=tensor_train)
        lossLabel2_val = sess.run(lossLabel2, feed_dict=tensor_train)
        lossMMD_val = sess.run(lossMMD, feed_dict=tensor_train)
        
        print(f"{i:<10} {loss_val:<15.6f} {lossLabel2_val:<15.6f} {lossMMD_val:<15.6f}")
        
        if i < train_steps:
            tensor_train = prepare_training_batch()



