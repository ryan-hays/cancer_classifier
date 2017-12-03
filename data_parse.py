"""
File to parse studies on cBioPortal.
First, we want to extract the most mutated genes in the studies by their mutations/nucleotide values.
"""
import sys
import os 
import numpy as np

def parse_mutation_file(file, gene_dict):
    t_barcode_index = file[0].index('Tumor_Sample_Barcode')
    visited_samples = set()

    for i in file:
        if i[0] not in gene_dict:
            gene_dict[i[0]] = 1
            visited_samples.add(i[t_barcode_index])
        else:
            if i[t_barcode_index] not in visited_samples:
                gene_dict[i[0]] += 1
    return gene_dict

def get_genes_for_cancer_type(files, k=10):
    """
    Parses each data file and returns k-most prevelant genes
    """
    gene_dict = {}
    for data in files:
        print(data)
        try:
            open_file = readdata(data)
            gene_dict = parse_mutation_file(open_file, gene_dict=gene_dict)
        except:
            print("Study not parsed")
            pass
    #Sort and extract top k
    sorted_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse = True)[:k]
    gene_list = [gene[0] for gene in sorted_genes]
    return gene_list


def get_cancer_data_files(mutation_data_name):
    """
    Gets a dictionary mapping cancer type to a list of data files

    Parameters: 
    mutation_data_name: (str) the file name for the mutation data. Consistent across all
    cancer types

    Returns:
    cancer_types_dict: (dict) cancer : list of data files
    """
    cancer_types_dict = {} #key = cancer_type, value = mutation data files
    cancer_types_folder = os.listdir("../cancer_classifier_data")
    #Iterate through every cancer type
    for cancer_type in cancer_types_folder:
        cancer_types_dict[cancer_type] = []
        for root, dirs, files in os.walk(os.path.join("../cancer_classifier_data", cancer_type)):
            if mutation_data_name in files:
                cancer_types_dict[cancer_type].append(os.path.join(root, mutation_data_name)) #collect all files of data for each cancer type
    return cancer_types_dict

def readdata(file):
    data = []
    for line in open(file,'r'):
        line = line.rstrip()
        info = line.split("\t")
        data.append(info)
    
    start_i = None

    for i in range(len(data)):
        if data[i][0] == "Hugo_Symbol":
            start_i = i
            break
    if start_i != None:
        data = data[start_i:]

    return data

def gene_features():
    """
    Select genes to be used as features

    Returns:
    gene_set: (list) contains HUGO names of the genes being used as the features 
    """
    cancer_type_dict = get_cancer_data_files(mutation_data_name="data_mutations_extended.txt")
    gene_set = []
    #Pull genes for each cancer type
    for cancer in cancer_type_dict:
        cancer_files = cancer_type_dict[cancer]
        genes = get_genes_for_cancer_type(cancer_files)
        for gene in genes:
            if gene not in gene_set: #make sure gene hasn't already been selected in another type
                gene_set.append(gene)
    return gene_set

def get_cancer_types(file_paths_by_cancer):
    """
    Maps cancer type to an index value 

    Parameters:
    file_paths_by_cancer: (dict) cancer : list of data files 

    Returns:
    cancer_dict: (dict) cancer : integer identifier
    """
    cancer_dict = {}
    cancer_index = 0
    for cancer in file_paths_by_cancer:
        cancer_dict[cancer] = cancer_index
        cancer_index += 1
    return cancer_dict

def get_samples(file_paths_by_cancer, cancer_dict, feature_vector_length):
    """
    Intializes a feature vector for samples in dataset and maps samples to their cancer type 

    Parameters:
    file_paths_by_cancer: (dict) cancer : list of data files 
    cancer_dict: (dict) cancer : integer identifier
    feature_vector_length: (int) length of feature vector

    Returns:
    training_set: (dict) sample : initialized feature vector 
    training_set_targets: (dict) sample : cancer integer identifer
    """
    training_set = {}
    training_set_targets = {}
    #Intialize training set and training set targets with all zeros
    for cancer in file_paths_by_cancer:
        list_of_files = file_paths_by_cancer[cancer]
        for f in list_of_files:
            try:
                data = readdata(f)
                t_barcode_index = data[0].index('Tumor_Sample_Barcode')
                for i in data[1:]:
                    if i[t_barcode_index] not in training_set: # and i[t_barcode_index] != 'Tumor_Sample_Barcode':
                        sample_id = i[t_barcode_index]
                        training_set[sample_id] = [0 for i in range(feature_vector_length)]
                        training_set_targets[sample_id] = cancer_dict[cancer]
            except Exception as e:
                print(e)
                pass
    return training_set, training_set_targets

def get_genes_for_samples_binary(gene_set, training_set, file_paths_by_cancer):
    """
    Builds a binary feature vector for each sample in the data set. Gene is "1" if mutated in
    given sample; "0" else

    Parameters:
    gene_set: (list) contains HUGO names of the genes being used as the features 
    training_set: (dict) sample : feature vector, where "feature vector" contains mutation data
    file_paths_by_cancer: (dict) cancer : list of data files 

    Returns:
    training_set: mutated "training_set" parameters where genes are set to "1" if mutated in sample
    """
    for cancer in file_paths_by_cancer:
        list_of_files = file_paths_by_cancer[cancer]
        for f in list_of_files:
            try:
                data = readdata(f)
                t_barcode_index = data[0].index('Tumor_Sample_Barcode')
                for i in data[1:]:
                    sample_id = i[t_barcode_index]
                    gene = i[0]
                    if gene in gene_set:
                        training_set[sample_id][gene_set.index(gene)] = 1
            except:
                #TODO fix error handling
                pass
    return training_set

def save_data_set(training_set, training_set_targets, dataset_name):
    """
    Saves training set feature vectors and target values as numpy array

    Parameters:
    training_set: (dict) sample : feature vector, where "feature vector" contains mutation data
    training_set_targets: (dict) sample : cancer integer identifer
    dataset_name: (str) dataset name for saving numpy array
    """
    sample_list = training_set.keys()
    training_set_final = [training_set[val] for val in sample_list]
    training_set_targets_final = [training_set_targets[val] for val in sample_list]
    np.save(dataset_name+".npy", np.array(training_set_final))
    np.save(dataset_name+"_targets"+".npy", np.array(training_set_targets_final))

def feature_build(gene_set):
    file_paths_by_cancer = get_cancer_data_files(mutation_data_name="data_mutations_extended.txt")
    cancers = get_cancer_types(file_paths_by_cancer) #Find cancers
    samples, training_set_targets = get_samples(file_paths_by_cancer, 
                                                cancer_dict=cancers, 
                                                feature_vector_length=len(gene_set)) #Build list of samples    
    training_set = get_genes_for_samples_binary(gene_set, samples, file_paths_by_cancer) #Update feature vectors 
    print(training_set)
    #Save targets
    save_data_set(training_set, training_set_targets, dataset_name="TEST")

def main():
    genes = gene_features() #select genes that we will use as gene features 
    # feature_build(gene_set=genes)#Now build feature vector for each sample in study 

if __name__ == "__main__":
    main()



