from bioinfokit import analys, visuz
import pandas as pd
import numpy as np
import json
import os
import pickle
import logging
from genedrugrecommender.utils import get_logger
from genedrugrecommender.disgenet import DisGenet
from genedrugrecommender.druggeneinteraction import DrugGeneInteraction
from genedrugrecommender.lincsregulatedb import LincsRegulateDb
from genedrugrecommender.limma import Limma
from genedrugrecommender.afimap import AfiMap
from genedrugrecommender.dataimport import DataImport

class DrugRecommenderConfig:
    def __init__(self,
                 data_folder, 
                 num_disease_suggestions,
                 num_drugs_suggestions,
                 num_top_genes):
        self.data_folder = data_folder
        self.num_disease_suggestions = num_disease_suggestions
        self.num_drugs_suggestions = num_drugs_suggestions
        self.num_top_genes = num_top_genes

class DrugRecommender:
    def __init__(self, config_filename):
        with open(config_filename, ) as f:
            config_json = json.load(f)
        self.logger = get_logger('DrugRecommender')
        self.config = DrugRecommenderConfig(**config_json)
        disgenet_db_path = os.path.join(self.config.data_folder, "disgenet_2020.db")
        self.disgenet = DisGenet(disgenet_db_path)
        lfile = os.path.join(self.config.data_folder, "LINCS_L1000")
        self.lincs = LincsRegulateDb(lfile)
        self.method = Limma()
        self.data_source = DataImport()

    def _get_best_drug_and_remove(self, all_drugs, gene_drugs):
        drug_scores = {}
        for drug in all_drugs:
            if drug not in drug_scores:
                drug_scores[drug] = 0
            for gene_set in gene_drugs:
                if drug in gene_set:
                    drug_scores[drug] += 1
        mx = -1
        drug = ''
        for drug_name, score in drug_scores.items():
            if mx < score:
                mx = score
                drug = drug_name
        for gene_set in gene_drugs:
            if drug in gene_set: 
                gene_set.remove(drug)
        return drug, mx
    
    def split_genes(self, top_genes, logFC):
        gene_drugs = []
        all_drugs = []

        for idx, gene in enumerate(top_genes):
            if logFC[idx] >= 0:
                downs = self.lincs.get_perturb_drugs(gene, "down")
                all_drugs.extend(downs)
                gene_drugs.append(set(downs))
            else:
                ups = self.lincs.get_perturb_drugs(gene, "up") 
                all_drugs.extend(ups)
                gene_drugs.append(set(ups))
        all_drugs = set(all_drugs)
        return all_drugs, gene_drugs

    def run(self, experiment_name):
        control, perturbed, genes = self.data_source.get_data_and_filter(experiment_name)
        # get the differentially expressed genes
        genes, p_values, logFC, sort_values = self.method.run(control, perturbed, genes)
        df = pd.DataFrame({'log2FC': logFC, 'p-value': p_values, 'GeneNames': genes})
        top_genes = genes[:self.config.num_top_genes]
        self.logger.info(f"=== The top {self.config.num_top_genes} most significant genes ===")
        for idx, gene in enumerate(top_genes):
            self.logger.info(f"{idx + 1} Gene: {gene} score {sort_values[idx]}")
        self.logger.info("====================================================")

        self.logger.info("doing volcano plot for show")
        
        visuz.gene_exp.volcano(df=df, 
                               geneid="GeneNames",
                               genenames = tuple(top_genes),
                               lfc='log2FC',
                               pv='p-value',
                               dim=(10, 10),
                               show=True,
                               gstyle=2,
                               figname=f"{experiment_name}_volcano")
        
        self.logger.info("saving volcano plot")
        
        visuz.gene_exp.volcano(df=df, 
                               geneid="GeneNames",
                               genenames = tuple(top_genes),
                               lfc='log2FC',
                               pv='p-value',
                               dim=(10, 10),
                               show=False,
                               gstyle=2,
                               figname=f"{experiment_name}_volcano")

        all_drugs, gene_drugs = self.split_genes(top_genes, logFC)
        
        self.logger.info("=== Top Suggested drugs ===")
        for idx in range(1, self.config.num_drugs_suggestions+1):
            drug_name, drug_score = self._get_best_drug_and_remove(all_drugs, gene_drugs)
            self.logger.info(f"{idx} Drug: {drug_name} score:{drug_score}")
        self.logger.info("===========================")
        self.logger.info("detecting possible diseases")
        # sugest actual disease
        sugestions = self.disgenet.get_possible_diseases_from_gene_list(top_genes)
        top_diseases = sugestions[:self.config.num_disease_suggestions]
        self.logger.info("=== The top 3 most-likely diseases ===")
        for sugestion in top_diseases:
            self.logger.info(f"{sugestion['disease_name']} score: {sugestion['disease_score']}")
        self.logger.info("======================================")
        self.logger.info("done")