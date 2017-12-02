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


def extract_data_files():
    """
    Returns a dictionary mapping cancer type to a list of data_mutation_files
    """
    name = "data_mutations_extended.txt"
    cancer_types_dict = {} #key = cancer_type, value = mutation data files
    cancer_types_folder = os.listdir("../cancer_classifier_data")
    #Iterate through every cancer type
    for cancer_type in cancer_types_folder:
        cancer_types_dict[cancer_type] = []
        for root, dirs, files in os.walk(os.path.join("../cancer_classifier_data", cancer_type)):
            if name in files:
                cancer_types_dict[cancer_type].append(os.path.join(root, name)) #collect all files of data for each cancer type
    return cancer_types_dict

def readdata(file):
    data = []
    for line in open(file,'r'):
        line = line.rstrip()
        info = line.split("\t")
        data.append(info)
    del data[0]
    return data

def gene_features():
    cancer_type_dict = extract_data_files()
    gene_features = []
    #Pull genes for each cancer type
    for cancer in cancer_type_dict:
        cancer_files = cancer_type_dict[cancer]
        genes = get_genes_for_cancer_type(cancer_files)
        for gene in genes:
            if gene not in gene_features: #make sure gene hasn't already been selected in another type
                gene_features.append(gene)
    return gene_features, cancer_type_dict

def feature_build(gene_set, file_paths_by_cancer):
    #Build list of samples
    cancer_dict = {} #maps cancer to an index
    cancer_index = 0
    training_set = {}
    training_set_targets = {}
    #Intialize training set and training set targets with all zeros
    for cancer in file_paths_by_cancer:
        cancer_dict[cancer] = cancer_index
        cancer_index += 1
        list_of_files = file_paths_by_cancer[cancer]
        for f in list_of_files:
            try:
                data = readdata(f)
                t_barcode_index = data[0].index('Tumor_Sample_Barcode')
                for i in data[1:]:
                    if i[t_barcode_index] not in training_set: # and i[t_barcode_index] != 'Tumor_Sample_Barcode':
                        sample_id = i[t_barcode_index]
                        training_set[sample_id] = [0 for i in range(len(gene_set))]
                        training_set_targets[sample_id] = cancer_dict[cancer]
            except Exception as e:
                print(e)
                pass
    #print("Here", training_set)
    #Update when finding genes
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
                pass
    #Concert to array 
    #print(len(training_set))
    print training_set
    print training_set_targets
    sample_list = training_set.keys()
    training_set_final = [training_set[val] for val in sample_list]
    training_set_targets_final = [training_set_targets[val] for val in sample_list]
    np.save('v2_training_set.npy', np.array(training_set_final))
    np.save('v2_training_set_targets.npy', np.array(training_set_targets_final))
    #print(gene_set)





if __name__ == "__main__":
    genes, cancer_types_dict = gene_features() #select genes that we will use as gene features 
    feature_build(gene_set=genes, file_paths_by_cancer=cancer_types_dict)#Now build feature vector for each sample in study 



