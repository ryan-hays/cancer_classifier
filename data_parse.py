"""
File to parse studies on cBioPortal.
First, we want to extract the most mutated genes in the studies by their mutations/nucleotide values.
"""
import sys


def parse_mutation_file(file):
    t_barcode_index = file[0].index('Tumor_Sample_Barcode')
    gene_dict = {}
    visited_samples = set()

    for i in file:
        if i[0] not in gene_dict:
            gene_dict[i[0]] = 1
            visited_samples.add(i[t_barcode_index])
        else:
            if i[t_barcode_index] not in visited_samples:
                gene_dict[i[0]] += 1


    sorted_genes = sorted(gene_dict.items(), key=lambda x: x[1], reverse = True)
    top = sorted_genes[:100]
    print sorted_genes
    print top
    






def readdata(file):
    data = []
    for line in open(file,'r'):
        line = line.rstrip()
        info = line.split("\t")
        data.append(info)
    del data[0]
    return data


def main():

    file_name = sys.argv[1]
    open_file = readdata(file_name)
    print parse_mutation_file(open_file)




if __name__ == "__main__":
    main()