"""
File to parse studies on cBioPortal.
First, we want to extract the most mutated genes in the studies by their mutations/nucleotide values.
"""
import sys
import os 

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

def get_genes_for_cancer_type(files, k=100):
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
    return sorted_genes


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

def main():
    cancer_type_dict = extract_data_files()
    gene_features = []
    #Pull genes for each cancer type
    for cancer in cancer_type_dict:
        cancer_files = cancer_type_dict[cancer]
        genes = get_genes_for_cancer_type(cancer_files)
        for gene in genes:
            if gene not in gene_features: #make sure gene hasn't already been selected in another type
                gene_features.append(gene)
    

# def get_web_data(cancer):
#     #cancer_studies #Say that this is our list of cancer studies
#     test = urllib2.urlopen("http://www.cbioportal.org/webservice.do?cmd=getCaseLists&cancer_study_id=%s") % (cancer)
#     print (test.read())

if __name__ == "__main__":
    main()