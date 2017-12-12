import tensorflow as tf  
import numpy as np 
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

tf.logging.set_verbosity(tf.logging.INFO)
# #Filepath for np arrays


def run_training(width, depth):
    training_set_file_features = "CNV_TEST_diff_916.npy"
    training_set_file_targets = "CNV_TEST_targets_diff_916.npy"
    K = 10

    raw_training_set_features = np.load(training_set_file_features)
    raw_training_set_targets = np.load(training_set_file_targets)

    list_training_set_features = np.array_split(raw_training_set_features, K)
    list_training_set_targets = np.array_split(raw_training_set_targets, K)

    print("SIZE", raw_training_set_features.shape[0]/10)
    accuracies = []
    hidden_unit_descript = []
    for i in range(1, depth+1):
        if i == 1:
            hidden_unit_descript.append(width)
        else:
            val = (3**i)*width/(5**i)
            if val == 0:
                val = 1
            hidden_unit_descript.append(val)
    print(hidden_unit_descript, "w")
    #K fold cross validation
    for i in range(len(list_training_set_features)):
        test_set_features = list_training_set_features[i]
        test_set_targets = list_training_set_targets[i]
        training_set_features = np.concatenate([list_training_set_features[j] for j in range(len(list_training_set_features)) if j != i])
        training_set_targets = np.concatenate([list_training_set_targets[j] for j in range(len(list_training_set_targets)) if j != i])
        print(training_set_features.shape)
        #Training set data and target values saved as numpy arrays
        NUM_OF_GENES = training_set_features.shape[1] #, or something like this 

        # Specify that all features have real-value data
        feature_columns = [tf.feature_column.numeric_column("genes", shape=[NUM_OF_GENES])] #num of genes

        # Build 3 layer DNN with 10, 20, 10 units respectively.
        classifier = tf.estimator.DNNClassifier(feature_columns=feature_columns,
                                                hidden_units=hidden_unit_descript,
                                                n_classes=9,
                                                dropout=None,
                                                model_dir="")

        # Define the training inputs
        train_input_fn = tf.estimator.inputs.numpy_input_fn(
            x={"genes": training_set_features},
            y=training_set_targets,
            num_epochs=10,
            batch_size=100,
            shuffle=True)

        train_input_fn_2 = tf.estimator.inputs.numpy_input_fn(
            x={"genes": training_set_features},
            y=training_set_targets,
            num_epochs=1,
            batch_size=100,
            shuffle=True)

        test_input_fn = tf.estimator.inputs.numpy_input_fn(
            x={"genes": test_set_features},
            y=test_set_targets,
            num_epochs=1,
            shuffle=False)

        # summary_hook = tf.train.SummarySaverHook(
        #     100,
        #     output_dir='training_summaries/',
        #     summary_op=tf.train.Scaffold(summary_op=tf.summary.merge_all()))

        classifier.train(input_fn=train_input_fn, steps=200000), #hooks=[summary_hook])
        ev = classifier.evaluate(input_fn=test_input_fn)

        # Evaluate accuracy.
        accuracy_score = classifier.evaluate(input_fn=test_input_fn)["accuracy"]
        accuracies.append(accuracy_score)
        train_acc = classifier.evaluate(input_fn=train_input_fn_2)["accuracy"]

        tf.summary.scalar('accuracy', accuracy_score)
        tf.summary.scalar('train_acc', train_acc)

        # sess = tf.Session()
        # with tf.Session() as sess:
        #     # Merge all the summaries and write them out to /tmp/mnist_logs (by default)
        #     merged = tf.summary.merge_all()
        #     train_writer = tf.summary.FileWriter('training_summaries',
        #                                           sess.graph)


        print("\nTest Accuracy: {0:f}\n".format(accuracy_score))
        print("\nTrain Accuracy: {0:f}\n".format(train_acc))

    print("\nAvg Accuracy: {0:f}\n".format(float(sum(accuracies))/float(len(accuracies))))
    avg_accuracy = float(sum(accuracies))/float(len(accuracies))
    return avg_accuracy

#Run for multiple widths/depth
if __name__ == "__main__":
    x_vals = []
    y_vals = []
    z_vals = []
    for width in [10, 50, 100, 250, 500, 1000]:
        for depth in range(1, 11):
            acc = run_training(width=width, depth=depth)
            #acc = 1
            x_vals.append(width)
            y_vals.append(depth)
            z_vals.append(acc)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(np.array(x_vals), np.array(y_vals), np.array(z_vals))
    #plt.plot_surface(x_vals, y_vals, z_vals)
    plt.show()









