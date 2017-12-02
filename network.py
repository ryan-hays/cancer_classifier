import tensorflow as tf  
import numpy as np 
import random

#Test this: create 100 samples, 10 features, 8 classes
# targets = []
# training_set = []
# for i in range(100):
# 	feature_vector = []
# 	for j in range(10):
# 		#Create feature 
# 		feature_vector.append(random.uniform(0.0, 1.0))
# 	training_set.append(feature_vector)
# 	targets.append(random.randint(0, 7))

# #Filepath for np arrays
training_set_file_features = "training_set_features.npy"
training_set_file_targets = "training_set_targets.npy"

# #Save test arrays
# np.save(training_set_file_features, np.array(training_set))
# np.save(training_set_file_targets, np.array(targets))


# #Training set data and target values saved as numpy arrays
training_set_features = np.load(training_set_file_features)
training_set_targets = np.load(training_set_file_targets)
NUM_OF_GENES = 10 #training_set_features.size, or something like this 

# Specify that all features have real-value data
feature_columns = [tf.feature_column.numeric_column("genes", shape=[NUM_OF_GENES])] #num of genes

# Build 3 layer DNN with 10, 20, 10 units respectively.
classifier = tf.estimator.DNNClassifier(feature_columns=feature_columns,
                                        hidden_units=[10, 20, 10],
                                        n_classes=8,
                                        model_dir="")

# Define the training inputs
train_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"genes": training_set_features},
    y=training_set_targets,
    num_epochs=None,
    shuffle=True)

test_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"genes": training_set_features},
    y=training_set_targets,
    num_epochs=1,
    shuffle=True)

classifier.train(input_fn=train_input_fn, steps=2000)
ev = classifier.evaluate(input_fn=test_input_fn)

# Evaluate accuracy.
accuracy_score = classifier.evaluate(input_fn=test_input_fn)["accuracy"]

print("\nTest Accuracy: {0:f}\n".format(accuracy_score))