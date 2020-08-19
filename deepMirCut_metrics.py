import tensorflow as tf
from keras import backend as K

def avg_proximity_metric():

    def prox(y_true, y_pred):
        eps = K.epsilon()
        beta = 1e10
        #print("y_pred SHAPE=")
        #print(K.int_shape(y_pred));
        #y_pred = tf.Print(y_pred,[y_pred], "y_pred=")
        #print("y_true SHAPE=")
        #print(K.int_shape(y_true));
        #y_true = tf.Print(y_true,[y_true], "y_true=")
        
        y_pred_T = tf.transpose(y_pred,perm=[0,2,1])
        y_true_T = tf.transpose(y_true,perm=[0,2,1])
        #print("y_pred_T SHAPE=")
        #print(K.int_shape(y_pred_T));
        #print("y_true_T SHAPE=")
        #print(K.int_shape(y_true_T));

        
        #y_range = tf.range(y_pred_T.shape.as_list()[-1], dtype=y_pred_T.dtype)
        y_pred_argmax = K.cast(K.argmax(y_pred_T,axis=-1),dtype='float32')
        y_true_argmax = K.cast(K.argmax(y_true_T,axis=-1),dtype='float32')
        #print("y_pred_argmax SHAPE=")
        #print(K.int_shape(y_pred_argmax))
        #print("y_true_argmax SHAPE=")
        #print(K.int_shape(y_true_argmax))        
        #y_pred_argmax = tf.Print(y_pred_argmax,[y_pred_argmax,y_pred_argmax[0].shape], "y_pred_argmax=")
        #y_true_argmax = tf.Print(y_true_argmax,[y_true_argmax], "y_true_argmax=")

        y_pred_slice = y_pred_argmax[:,1:5]
        y_true_slice = y_true_argmax[:,1:5]
        #y_pred_argmax = tf.Print(y_pred_argmax,[y_pred_argmax,y_pred_argmax[0].shape,y_pred_slice,y_pred_slice[0].shape], "y_pred_argmax_slice=")

        y_prox = K.abs(y_pred_slice - y_true_slice)
        #y_prox = tf.Print(y_prox,[y_prox],"y_prox=");
        y_prox_sum = K.sum(y_prox,axis=-1)
        
        return K.mean(y_prox_sum, axis=-1)
    
    return prox
