"""
File to parse studies on cBioPortal.
First, we want to extract the most mutated genes in the studies by their mutations/nucleotide values.
"""
import sys
import os 
import numpy as np

def CNV_parse(file, gene_dict):

    for line in file[1:]:
        if line[0] in gene_dict:
            gene_dict[line[0]] += sum([(i**2)**(1/2) for i in line[1:]])
        else:
            gene_dict[line[0]] = sum([(i**2)**(1/2) for i in line[1:]])
    return gene_dict

def get_genes_for_cancer_type(files, k=100):
    """
    Parses each data file and returns k-most prevelant genes
    """
    gene_dict = {}
    for data in files:
        print(data)
        try:
            open_file = readdata(data, True)
            gene_dict = CNV_parse(open_file, gene_dict=gene_dict)
        except:
            print "Study not parsed"

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
    cancer_types_dict.pop(".DS_Store", None)
    return cancer_types_dict

def readdata(file, CNV = False):
    data = []
    for line in open(file,'r'):
        line = line.rstrip()
        info = line.split("\t")
        data.append(info)
    
    start_i = None

    #Removes text above header if there is any
    for i in range(len(data)):
        if data[i][0] == "Hugo_Symbol":
            start_i = i
            break
    if start_i != None:
        data = data[start_i:]

    if not CNV:
        return data

    else:
        #Removes useless information, making the data easier to use
        data = remove_index(data, "Entrez_Gene_ID")
        data = remove_index(data, "Entrez_Gene_Id")
        data = remove_index(data, "Cytoband")

        #Removes genes with incomplete data
        new_data = []
        new_data.append(data[0])
        for line in data[1:]:
            try:
                symbol = line[0]
                nums = [int(i) for i in line[1:]]
                new_data.append([symbol]+nums)
            except:
                pass

        return new_data



def remove_index(data, value):
    """
    Removes an index in the data for every gene
    """
    in_data = False
    val_index = None
    if value in data[0]:
        in_data = True
        val_index = data[0].index(value)

    if in_data:
        for line in data:
            del line[val_index]

    return data 



def gene_features():
    """
    Select genes to be used as features

    Returns:
    gene_set: (list) contains HUGO names of the genes being used as the features 
    """
    cancer_type_dict = get_cancer_data_files(mutation_data_name="data_CNA.txt")
    gene_set = []
    #Pull genes for each cancer type
    for cancer in cancer_type_dict:
        cancer_files = cancer_type_dict[cancer]
        if len(cancer_files) == 0:
            continue
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
    print cancer_dict
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
                data = readdata(f, True)
                
                for i in data[0][1:]: 
                    if i not in training_set: 
                        sample_id = i
                        training_set[sample_id] = [0 for i in range(feature_vector_length)]
                        training_set_targets[sample_id] = cancer_dict[cancer]
            except Exception as e:
                print(e)
                pass
    return training_set, training_set_targets

def get_CNV_for_samples(gene_set, training_set, file_paths_by_cancer):
    """
    Builds feature vector for each sample in the data set. Location is the CNV.

    Parameters:
    gene_set: (list) contains HUGO names of the genes being used as the features 
    training_set: (dict) sample : feature vector, where "feature vector" contains mutation data
    file_paths_by_cancer: (dict) cancer : list of data files 

    Returns:
    training_set: CNV "training_set" parameters where genes are set to CNV value 
    """
    gene_fast = set(gene_set)
    for cancer in file_paths_by_cancer:
        list_of_files = file_paths_by_cancer[cancer]
        for f in list_of_files: #For each file 
            print f

            data = readdata(f, True)
            
            for i in data[1:]: #For each gene in file 

                gene = i[0]
                if gene in gene_fast: #if gene is one of the top genes
                    for j in range(1, len(i[1:])): #Assign CNV value for of sample for a gene
                        sample_id = data[0][j]
                        training_set[sample_id][gene_set.index(gene)] = i[j]

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
    file_paths_by_cancer = get_cancer_data_files(mutation_data_name="data_CNA.txt")
    cancers = get_cancer_types(file_paths_by_cancer) #Find cancers
    samples, training_set_targets = get_samples(file_paths_by_cancer, 
                                                cancer_dict=cancers, 
                                                feature_vector_length=len(gene_set)) #Build list of samples 
    training_set = get_CNV_for_samples(gene_set, samples, file_paths_by_cancer) #Update feature vectors 
    print training_set
    #Save targets
    save_data_set(training_set, training_set_targets, dataset_name="TEST")

def main():
    genes = gene_features() #select genes that we will use as gene features
    # genes = ['PKIG', 'SERINC3', 'TTPAL', 'ADA', 'SEMG2', 'MATN4', 'HM13', 'TM9SF4', 'FITM2', 'TPX2', 'R3HDML', 'HNF4A', 'SLPI', 'RBPJL', 'FOXS1', 'MYLK2', 'YWHAB', 'PI3', 'HCK', 'WFDC5', 'RIMS4', 'KCNS1', 'SEMG1', 'JPH2', 'SYS1', 'ASIP', 'PABPC1L', 'KCNK15', 'GDAP1L1', 'POFUT1', 'DBNDD2', 'WFDC12', 'SDC4', 'REM1', 'TP53TG5', 'TTLL9', 'DEFB124', 'STK4', 'DUSP15', 'PDRG1', 'WISP2', 'PLAGL2', 'COX4I2', 'ID1', 'XKR7', 'KIF3B', 'ASXL1', 'ITCH', 'EIF2S2', 'PIGT', 'TOMM34', 'DEFB123', 'DEFB122', 'BCL2L1', 'WFDC2', 'AHCY', 'TOX2', 'PIGU', 'MAP1LC3A', 'DYNLRB1', 'NCOA6', 'WFDC9', 'DEFB118', 'DEFB119', 'WFDC10B', 'WFDC10A', 'RALY', 'WFDC11', 'WFDC13', 'TP53INP2', 'DEFB121', 'SPINT3', 'C20orf112', 'MYH7B', 'TRPC4AP', 'PROCR', 'WFDC8', 'WFDC6', 'NFS1', 'GTSF1L', 'SPINT4', 'RBM39', 'ROMO1', 'EDEM2', 'GGT7', 'ERGIC3', 'CHMP4B', 'ACSS2', 'WFDC3', 'CEP250', 'GSS', 'C20orf173', 'IFT52', 'MMP24', 'DNTTIP1', 'UBE2C', 'MYBL2', 'ACOT8', 'SNX21', 'TNNC2', 'CDKN2A', 'CDKN2B', 'C9orf53', 'MTAP', 'EGFR', 'SNORA73|ENSG00000252054.1', 'SEC61G', 'LANCL2', 'RN7SL151P', 'MIR31HG', 'VOPP1', 'SNORD39|ENSG00000264379.1', 'IFNE', 'MIR31', 'DMRTA1', 'IFNA1', 'VSTM2A', 'IFNA8', 'IFNA2', 'IFNA13', 'IFNA6', 'KLHL9', 'IFNA5', 'PAOX', 'FKBP9L', 'SPRN', 'MTG1', 'FRG2B', 'SYCE1', 'CYP2E1', 'MIR3944', 'ECHS1', 'MIR202', 'CALY', 'FUOM', 'ZNF511', 'ADAM8', 'PRAP1', 'VENTX', 'KNDC1', 'UTF1', 'TUBGCP2', 'INPP5A', 'GPR123', 'IFNA14', 'TTC40', 'PTEN', 'JAKMIP3', 'STK32C', '14-Sep', 'IFNA17', 'EBF3', 'TCERG1L', 'MIR378C', 'LRRC27', 'CTBP2', 'PWWP2B', 'C10orf91', 'snoU13|ENSG00000238354.1', 'GLRX3', 'DPYSL4', 'LINC00959', 'CPA2', 'CPA1', 'CPA4', 'CPA5', 'MIR4297', 'CEP41', 'MKLN1', 'PTPRE', 'snoU13|ENSG00000238336.1', 'snoU13|ENSG00000239044.1', 'RNA5SP246', 'MIR29A', 'MIR29B1', 'CLRN3', 'COPG2', 'C10orf90', 'MGMT', 'TSGA13', 'IFNA4', 'KLF14', 'FAM196A', 'TMEM209', 'MIR335', 'TEX36', 'MEST', 'METTL10', 'NPS', 'MKI67', 'SSMEM1', 'PPP2R2D', 'TIAL1', 'FAM53B', 'INPP5F', 'RN7SL846P', 'TACC2', 'DOCK1', 'FOXI2', 'MIR4682', 'TERT', 'NKD2', 'CLPTM1L', 'SLC6A3', 'C5orf38', 'IRX2', 'IRX1', 'DNAH5', 'TRIO', 'FAM105A', 'FAM105B', 'ANKH', 'C5orf49', 'MTRR', 'FASTKD3', 'ADAMTS16', 'SEMA5A', 'LOC285692', 'NSUN2', 'TAS2R1', 'SRD5A1', 'ADCY2', 'MED10', 'FAM173B', 'UBE2QL1', 'MYO10', 'SNORD123', 'FLJ33360', 'CCT5', 'CMBL', 'KIAA0947', 'MARCH11', 'ZNF622', 'FAM134B', 'ROPN1L', 'CTNND2', 'MARCH6', 'BASP1', 'LOC285696', 'DAP', 'CDH18', 'PMCHL1', 'CDH12', 'PRDM9', 'RANBP3L', 'SKP2', 'LMBRD2', 'C5orf22', 'CDH10', 'NIPBL', 'UGT3A2', 'IL7R', 'C5orf42', 'CDH6', 'SLC1A3', 'CAPSL', 'RICTOR', 'UGT3A1', 'ZFR', 'CDH9', 'SPEF2', 'NPR3', 'PDZD2', 'FYB', 'OSMR', 'SUB1', 'MTMR12', 'DAB2', 'GOLPH3', 'C9', 'CARD6', 'NUP155', 'WDR70', 'TTC33', 'ADAMTS12', 'PTGER4', 'EGFLAM', 'SNORD72', 'PRKAA1', 'PRLR', 'RPL37', 'RXFP3', 'TARS', 'HEATR7B2', 'GDNF', 'LIFR', 'TTC23L', 'AGXT2', 'C7', 'C6', 'RAI14', 'ZBTB7B', 'ADAM15', 'TDRD10', 'DNAJC21', 'RAD1', 'DCST2', 'DCST1', 'CNTNAP2', 'TMEM178B', 'GALNTL5', 'GBX1', 'PAXIP1', 'TMUB1', 'ACTR3C', 'AGAP3', 'IQCA1P1', 'SLC4A2', 'ABCB8', 'FASTK', 'BLACE', 'UBE3C', 'ESYT2', 'CNPY1', 'NOS3', 'GALNT11', 'KCNH2', 'DNAJB6', 'RARRES2', 'ASIC3', 'AOC1', 'VIPR2', 'SSPO', 'MNX1', 'EN2', 'RBM33', 'RNF32', 'REPIN1', 'CRYGN', 'LRRC61', 'MIR671', 'SHH', 'MIR3907', 'LMBR1', 'ZBED6CL', 'GIMAP8', 'GIMAP4', 'GIMAP5', 'GIMAP6', 'GIMAP7', 'GIMAP1', 'GIMAP2', 'ASB10', 'ABCF2', 'LINC00689', 'NCAPG2', 'SMARCD3', 'WDR60', 'ZNF467', 'ATG9B', 'TMEM176A', 'TMEM176B', 'HTR5A', 'WDR86', 'MIR595', 'RHEB', 'ATP6V0E2', 'PRKAG2', 'C7orf13', 'CHPF2', 'CDK5', 'ZNF775', 'ZNF862', 'PTPRN2', 'NUB1', 'INSIG1', 'MIR5707', 'DPP6', 'NOM1', 'ZNF212', 'CTAGE8', 'NOBOX', 'RNA5SP249', 'ARHGEF35', 'MIR548F4', 'ZNF425', 'ZNF783', 'ZNF786', 'CTAGE4', 'OR2A20P', 'TPK1', 'KRBA1', 'RNY5', 'RNY4', 'RNY1', 'EXOC4', 'ZNF767', 'EZH2', 'PLXNA4', 'RNY3', 'C7orf33', 'ZNF398', 'OR2A42', 'PDIA4', 'ARHGEF5', 'ZNF746', 'ZNF282', 'OR2A1', 'DUB3', 'NME1-NME2', 'LOC654342', 'LOC653712', 'LOC654346', 'LOC654349', 'LOC654350', 'LOC654348', 'LOC654347', 'LOC654345', 'LOC654341', 'LOC654343', 'LOC654340', 'LOC654335', 'LOC653786', 'LOC653513', 'LOC653653', 'LOC653827', 'LOC653822', 'LOC653814', 'LOC653501', 'LOC653820', 'LOC653816', 'LOC653809', 'LOC653818', 'LOC653813', 'LOC653815', 'LOC653812', 'LOC653810', 'LOC653811', 'LOC653803', 'LOC653807', 'LOC653801', 'LOC653826', 'LOC653825', 'LOC653806', 'LOC653808', 'LOC653805', 'LOC653824', 'LOC653823', 'LOC653802', 'LOC653804', 'LOC653799', 'LOC653800', 'LOC653798', 'LOC653796', 'LOC653795', 'LOC653821', 'LOC653817', 'LOC653791', 'LOC653793', 'LOC653797', 'LOC653792', 'LOC653789', 'LOC653794', 'LOC653790', 'LOC653788', 'LOC653782', 'LOC653783', 'LOC653779', 'LOC653780', 'LOC653784', 'LOC653771', 'LOC653770', 'LOC653776', 'LOC653769', 'LOC653777', 'LOC653768', 'LOC653781', 'LOC653773', 'LOC653772', 'LOC653774', 'LOC653760', 'LOC653764', 'LOC653762', 'LOC653775', 'LOC653767', 'LOC653766', 'LOC653763', 'LOC653761', 'LOC653759', 'LOC653746', 'LOC653758', 'LOC653757', 'LOC653745', 'LOC653752', 'LOC653787', 'LOC653751', 'LOC653753', 'LOC653736', 'LOC653750', 'LOC653748', 'LOC653749', 'LOC653747', 'LOC653740', 'LOC653742', 'LOC653739', 'LOC653743', 'LOC653741', 'LOC653744', 'MIR548D2', 'C3orf74', 'CCDC72', 'hsa-mir-566', 'C3orf45', 'SNORD19B', 'C3orf77', 'SNORD19', 'VENTXP7', 'LOC401052', 'hsa-mir-191', 'hsa-mir-1226', 'C3orf75', 'C3orf71', 'hsa-mir-425', 'C3orf24', 'TMEM111', 'GALNTL2', 'hsa-mir-563', 'hsa-mir-564', 'C3orf39', 'ZNF167', 'C3orf19', 'hsa-mir-885', 'ARPP-21', 'CCRL2', 'hsa-mir-128-2', 'MIR128-2', 'C3orf23', 'LOC646498', 'KBTBD5', 'SNORA62', 'C3orf32', 'hsa-let-7g', 'C3orf10', 'C3orf42', 'C3orf48', 'hsa-mir-26a-1', 'WDR51A', 'LAMB2L', 'C3orf54', 'CDH29', 'C3orf31', 'hsa-mir-138-1', 'DVWA', 'CCBP2', 'hsa-mir-135a-1', 'MIR138-1', 'LOH3CR2A', 'TUSC4', 'C3orf51', 'HESRG', 'C3orf63', 'FAM116A', 'PBRM1', 'RAD54L2', 'TEX264', 'VPRBP', 'GLT8D1', 'GPR62', 'SNORD69', 'ALAS1', 'PARP3', 'RBM15B', 'MANF', 'ACY1', 'RRP9', 'GNL3', 'PCBP4', 'ABHD14A', 'ABHD14B', 'RPL29', 'DUSP7', 'DOCK3', 'BAP1', 'TWF2', 'DNAH1', 'NISCH', 'TLR9', 'SEMA3G', 'SPCS1', 'TNNC1', 'PHF7', 'DAG1', 'NAT6', 'SLC38A3', 'NT5DC2', 'RNF123', 'NEK4', 'STAB1', 'MST1', 'AMIGO3', 'GRM2', 'GMPPB', 'GNAT1', 'HYAL3', 'HYAL1', 'GNAI2', 'IFRD2', 'SEMA3B', 'PLEKHA6', 'PLXNA2', 'MDM4', 'PIK3C2B', 'LRRN2', 'PPP1R15B', 'NFASC', 'CDK18', 'MIR205HG', 'LEMD1', 'MFSD4', 'MIR135B', 'LEMD1-AS1', 'CNTN2', 'KISS1', 'GOLT1A', 'KLHDC8A', 'SRGAP2', 'ELK4', 'NUCKS1', 'REN', 'SOX13', 'CAMK1G', 'RBBP5', 'SLC45A3', 'ETNK2', 'NUAK2', 'DSTYK', 'PM20D1', 'SLC41A1', 'MIR4260', 'ESRRG', 'LAMB3', 'CD46', 'AVPR1B', 'CR1L', 'HSD11B1', 'MIR5191', 'G0S2', 'FAM72A', 'RASSF5', 'PRELP', 'OPTC', 'IKBKE', 'LINC00303', 'FMOD', 'MAPKAPK2', 'C4BPA', 'DYRK3', 'CR1', 'BTG2', 'CHIT1', 'PIGR', 'LINC01136', 'FCAMR', 'PFKFB2', 'EIF2D', 'C4BPB', 'YOD1', 'TRAF3IP3', 'ADORA1', 'MYBPH', 'CR2', 'CHI3L1', 'ATP2B4', 'IPO9-AS1', 'NAV1', 'MYOG', 'SYT2', 'PPFIA4', 'LAX1', 'TMCC2', 'TMEM81', 'ZC3H11A', 'KLHL12', 'ADIPOR1', 'RABIF', 'USH2A', 'TMEM183A', 'KDM5B', 'CYB5R1', 'SLC26A9', 'SIPA1L2', 'IRF6', 'SNRPE', 'PPP1R12B', 'PKP1', 'CTSE', 'CD34', 'C1orf186', 'IL24', 'PTPRVP', 'IL10', 'LGR6', 'IL20', 'IL19', 'DIEXF', 'C1orf116', 'PTPN7', 'UBE2T']
    feature_build(gene_set=genes)#Now build feature vector for each sample in study 

if __name__ == "__main__":
    main()






