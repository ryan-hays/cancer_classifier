from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import matplotlib.pyplot as plt

training_set_file_features = "PCA_SOMA_TEST.npy"
training_set_file_targets = "PCA_SOMA_TEST_targets.npy"

# #Training set data and target values saved as numpy arrays
training_set_features = np.load(training_set_file_features)
print(training_set_features.shape)
training_set_targets = np.load(training_set_file_targets)
# test_set_features = np.load(training_set_file_features)[:500]
# test_set_targets = np.load(training_set_file_targets)[:500]

def plot_pca(pc):
	for i in range(3):
		x_vals = [val[i] for val in pc]
		y_vals = [val[i+1] for val in pc]
		k_means = KMeans(n_clusters = 7)
		cluster_labels = k_means.fit_predict(pc)
		plt.figure(i)
		plt.scatter(x_vals,y_vals,c=cluster_labels)
		plt.title('K means (clusters = 7)')
		plt.xlabel('PC'+str(i+1))
		plt.ylabel('PC'+str(i+2)) 
		plt.show()
		#plt.savefig('pc'+str(i+1)+'_clusters.png')

max_acc = 0.0

#Kmeans
while(True):
	#clusters = AgglomerativeClustering(n_clusters=7)
	clusters = PCA(n_components=2)
	pc_vals = clusters.fit_transform(training_set_features).tolist()
	plot_pca(pc_vals)
	break

	kmeans = KMeans(n_clusters=7)
	kmeans_predict = kmeans.fit_predict(pc_vals).tolist()
	#kmeans_predict = clusters.fit_predict(training_set_features).tolist()


	#pca_struct = PCA(n_components=4)





	#Convert to list for convience
	training_targets = training_set_targets.tolist()
	#Determiine distance
	true_groupings = {}
	for i in range(len(training_set_targets)):
		if training_set_targets[i] not in true_groupings:
			true_groupings[training_set_targets[i]] = [i]
		else:
			true_groupings[training_set_targets[i]].append(i)

	predict_groupings = {}
	for i in range(len(kmeans_predict)):
		prediction = kmeans_predict[i]
		if prediction not in predict_groupings:
			predict_groupings[prediction] = [i]
		else:
			predict_groupings[prediction].append(i)

	#Label mapping maps predict label to true label integer
	label_map = {}

	#Now distance: calculated by the percent of matching items
	for predict_group in predict_groupings:
		#Find best label
		predict_group_cluster = predict_groupings[predict_group]
		max_matching = -1
		match_group = None
		for true_group in true_groupings:
			if true_group in label_map.values():
				continue
			matching = 0
			true_group_cluster = true_groupings[true_group]
			for item in predict_group_cluster:
				if item in true_group_cluster:
					matching += 1
			#print(predict_group, true_group, matching, max_matching)
			if matching > max_matching:
				match_group = true_group
				max_matching = matching
		label_map[predict_group] = match_group
	#print(label_map)

	#Create new prediction list
	new_predictions = []
	for predict in kmeans_predict:
		new_predictions.append(label_map[predict])

	#Now test accuracy
	correct = 0
	for i in range(len(new_predictions)):
		if new_predictions[i] == training_set_targets[i]:
			correct += 1
	accuracy = float(correct)/float(len(new_predictions))
	print(accuracy)
	if accuracy > max_acc:
		max_acc = accuracy
		print("Max Found ", max_acc)





