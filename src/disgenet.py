import sqlite3

class DisGenet:
    
    def __init__(self, sqlite_filename):
        self.sqlite_filename = sqlite_filename
        self.con = sqlite3.connect(self.sqlite_filename)
    
    def get_gene_diseases(self, gene_name):
        cursor = self.con.cursor()
        result = cursor.execute(f"SELECT diseaseAttributes.diseaseName,\
        geneDiseaseNetwork.score\
        FROM geneDiseaseNetwork, geneAttributes, diseaseAttributes \
        WHERE geneAttributes.geneName like '{gene_name}'\
        AND geneDiseaseNetwork.diseaseNID = diseaseAttributes.diseaseNID\
        AND geneDiseaseNetwork.geneNID = geneAttributes.geneNID")
        return_list = []
        for row in result:
            return_list.append({'disease_name': row[0], 'score': row[1]})
        return sorted(return_list, key = lambda i: i['score'], reverse=True)

    def get_possible_diseases_from_gene_list(self, gene_list):
        cursor = self.con.cursor()
        disease_scores = {}
        for gene in gene_list:
            diseases = self.get_gene_diseases(gene)
            for disease in diseases:
                disease_name = disease['disease_name']
                disease_score = disease['score']
                if disease_name in disease_scores:
                    disease_scores[disease_name] += disease_score
                else:
                    disease_scores[disease_name] = disease_score
        return_list = []
        for disease_name, score in disease_scores.items():
            return_list.append({'disease_name':disease_name, 'disease_score':score})
        return sorted(return_list, key = lambda i: i['disease_score'], reverse=True)

    def get_disease_genes(self, disease_name):
        cursor = self.con.cursor()
        result = cursor.execute(f"SELECT geneAttributes.geneName,\
        geneDiseaseNetwork.score\
        FROM geneDiseaseNetwork, geneAttributes, diseaseAttributes \
        WHERE diseaseAttributes.diseaseName like '{disease_name}'\
        AND geneDiseaseNetwork.diseaseNID = diseaseAttributes.diseaseNID\
        AND geneDiseaseNetwork.geneNID = geneAttributes.geneNID")
        return_list = []
        for row in result:
            return_list.append({'gene_name': row[0], 'score': row[1]})
        return sorted(return_list, key = lambda i: i['score'], reverse=True)

def test_get_genes_from_disease(sqlite_file = 'data/disgenet_2020.db', disease_name = 'Abscess'):
    dg = DisGenet(sqlite_file)
    disease_genes = dg.get_disease_genes(disease_name)
    print(disease_genes)
    print(f"found {len(disease_genes)} genes associated to this disease")

def test_get_diseases_from_gene(sqlite_file = 'data/disgenet_2020.db', gene_name = 'A2M'):
    dg = DisGenet(sqlite_file)
    genes_diseases = dg.get_gene_diseases(gene_name)
    print(genes_diseases)
    print(f"found {len(genes_diseases)} diseases associated to this gene")

def test_get_possible_diseases_from_gene_list(sqlite_file = 'data/disgenet_2020.db', gene_list = ['A2M','A2MP1']):
    dg = DisGenet(sqlite_file)
    genes_diseases = dg.get_possible_diseases_from_gene_list(gene_list)
    print(genes_diseases)

# test_get_genes_from_disease()
# test_get_diseases_from_gene()
test_get_possible_diseases_from_gene_list()