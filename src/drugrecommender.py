from disgenet import DisGenet
from druggeneinteraction import DrugGeneInteraction
from drugrecommender import DrugRecommender
from lincsregulatedb import LincsRegulateDb
import json
import os

class DrugRecommenderConfig:
    def __init__(self, data_folder, num_disease_suggestions):
        self.data_folder = data_folder
        self.num_disease_suggestions = num_disease_suggestions

class DrugRecommender:
    def __init__(self, config_filename):
        config_json = json.load(config_filename)
        self.config = DrugRecommenderConfig(**config_json)
        disgenet_db_path = os.path.join(self.config.data_folder, "disgenet_2020.db")
        self.disgenet = DisGenet(disgenet_db_path)
        self.lincs = LincsRegulateDb(self.config.data_folder)
    
    # use methods to get up or down regulated genes from experiment data
    def get_DEG(self, control, perturbed):
        return []

    # plot differentially expressed genes
    def plot_DEG(self, degs):
        pass

    def split_DEG(self, degs):
        return [],[]

    def run(self, control, perturbed, gene_list):
        # get the differentially expressed genes
        degs = self.get_DEG(control, perturbed)
        # plot differentially expressed genes
        self.plot_DEG(degs)
        # split DEG results in up or down
        up_regulated, down_regulated = self.split_DEG(degs)
        # sugest actual disease
        suggestion_input = up_regulated + down_regulated
        sugestion = self.disgenet.get_possible_diseases_from_gene_list(suggestion_input)
        print(f"Printing top {self.config.num_disease_suggestions}")
        print(sugestion[:self.config.num_disease_suggestions])