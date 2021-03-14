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