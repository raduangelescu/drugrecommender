import numpy as np
import GEOparse
import pandas as pd
import numpy as np
import json
import os
import pickle
from genedrugrecommender.afimap import AfiMap

class DataImport:
    def __init__(self):
        self.afimap = AfiMap('data')
        self.folder = 'data'
    
    def _quantile_normalize(self, df_input):
        df = df_input.copy()
        dic = {}
        for col in df:
            dic.update({col: sorted(df[col])})
        sorted_df = pd.DataFrame(dic)
        rank = sorted_df.mean(axis=1).tolist()
        for col in df:
            t = np.searchsorted(np.sort(df[col]), df[col])
            df[col] = [rank[i] for i in t]
        return df

    def _log_if_necessary(self, data):
        if np.max(data) - np.min(data) > 100:
            data = np.where(data == 0, 1, data)
            return np.log2(data)
        return data

    def _repair_nan_fast(self, nparray):
        POS_INF = 10000000
        NEG_INF = 0
        mean = np.nanmean(nparray, axis=1)
        mean = np.nan_to_num(mean, copy=True, nan=0,
                             posinf=POS_INF,
                             neginf=NEG_INF)
        inds = np.where(np.isnan(nparray))
        nparray[inds] = np.take(mean, inds[0])
        return nparray

    def _get_data_from_GEO(self, geo_id = 'GSE15852'):
        cache_path = os.path.join("data", "geo_cache")
        gse = GEOparse.get_GEO(geo=geo_id, destdir=cache_path)
        control = ['Normal']
        phenotype_data = gse.phenotype_data
        columns = phenotype_data.columns
        info_experiment_idx = columns.get_loc('title')
        gsm_ids_idx = columns.get_loc('geo_accession')
        gsm_type = list(phenotype_data.values[:, info_experiment_idx])
        gsm_ids = list(phenotype_data.values[:, gsm_ids_idx])
        control_gsms = []
        perturbation_gsms = []
        raw_control_data = []
        raw_perturbed_data = []
        for idx in range(0, len(gsm_type)):
            gsm_id = gsm_ids[idx]
            table = gse.gsms[gsm_id].table
            value_idx = table.columns.get_loc('VALUE')
            values = gse.gsms[gsm_id].table.values[:, value_idx].tolist()
            if "Normal" in gsm_type[idx]:
                control_gsms.append(gsm_id)
                raw_control_data.append(values)
            else:
                perturbation_gsms.append(gsm_id)
                raw_perturbed_data.append(values)

        if not control_gsms:
            self.logger('[no col]no control for {geo_id}')
            return None

        genes_afi = gse.gsms[control_gsms[0]].table.values[:, 0]
        genes = self.afimap.convert_afi_list_to_gene_symbols(genes_afi)
        np_control_raw = np.array(raw_control_data)
        np_perturbed_raw = np.array(raw_perturbed_data)

        control = self._repair_nan_fast(np_control_raw)
        perturbed = self._repair_nan_fast(np_perturbed_raw)
        return control, perturbed, genes

    def _fix_bad_data(self, control_array, perturbed_array, genes):
        max_replicates = 3
        if len(control_array) < len(genes):
            # some experiments need to be transposed
            control = np.array(control_array).T.tolist()
            perturbed = np.array(perturbed_array).T.tolist()
        else:
            control = control_array
            perturbed = perturbed_array
        idx = 0
        # for control pick the first replicates
        control = [x[:max_replicates] for x in control]
        # for perturbed pick the last replicates
        perturbed = [x[-max_replicates:] for x in perturbed]
        # the above pick was done to favorize timeseries experiments
        control = self._log_if_necessary(np.array(control))
        perturbed = self._log_if_necessary(np.array(perturbed))

        control = self._quantile_normalize(pd.DataFrame(control))
        perturbed = self._quantile_normalize(pd.DataFrame(perturbed))

        return control, perturbed, genes

    def get_data_and_filter(self, experiment_name):
        file_name = f"{experiment_name}.pkl"
        path = os.path.join(self.folder, file_name)
        if os.path.isfile(path):
            with open(path, 'rb') as handle:
                print("loading data from cache")
                data = pickle.load(handle)
                control = data['control']
                perturbed = data['perturbed']
                genes = data['genes']
                return control, perturbed, genes

        print("getting data from GEO")
        control, perturbed, genes = self._get_data_from_GEO(experiment_name)
        print("fixing bad data")
        control, perturbed, genes = self._fix_bad_data(control, perturbed, genes)
        print("saving data")
        data_save = {
            'control':control,
            'perturbed':perturbed,
            'genes':genes
            }
        with open(path, 'wb') as handle:
            pickle.dump(data_save, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return control, perturbed, genes