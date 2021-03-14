import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import (pandas2ri, numpy2ri)
from rpy2.robjects.conversion import localconverter
import os
import numpy as np

class Limma:
    def __init__(self):
        current_file_path = os.path.dirname(os.path.abspath(__file__))
        self.abs_r_file_path = os.path.join(current_file_path,
                                            'limma.r')

    def sort(self, genes, values, logFC):
        genes = np.array(list(genes))
        abs_vals = np.absolute(values)
        max_value = np.max(abs_vals)
        sort_values = (np.ones(values.shape)* max_value - abs_vals)
        sort_values = -np.multiply(sort_values, np.absolute(logFC))
        sort_indices = np.argsort(sort_values)
        values = values[sort_indices]
        logFC = logFC[sort_indices]
        genes = genes[sort_indices]
        genes = [str(x) for x in genes]
        return genes, values, logFC, sort_values
    
    def prepare_data(self, A, B, genes):
        pdA = pd.DataFrame(A)
        pdB = pd.DataFrame(B)
        data_df = pd.concat([pdA, pdB],
                            axis=1,
                            ignore_index=True)
        mask = ([0] * pdA.shape[1])
        mask.extend([1] * pdB.shape[1])
        mask = pd.DataFrame(mask)
        return [data_df, mask, genes]

    def run(self, A, B, genes):
        pdA = pd.DataFrame(A)
        pdB = pd.DataFrame(B)
        data_df = pd.concat([pdA, pdB],
                            axis=1,
                            ignore_index=True)
        mask = ([0] * pdA.shape[1])
        mask.extend([1] * pdB.shape[1])
        mask = pd.DataFrame(mask)
        
        numpy2ri.activate()
        pandas2ri.activate()
        r = ro.r
        r.source(self.abs_r_file_path)
        data_df, mask, gene_df = self.prepare_data(A, B, genes)
        converters = ro.default_converter + pandas2ri.converter
        with localconverter(converters):
            r_data = ro.conversion.py2rpy(data_df)
            r_class = ro.conversion.py2rpy(mask)
            result = r.calculate_limma(r_data, r_class)
            result_py = ro.conversion.rpy2py(result)
            p_values = np.array(list(result_py['adj.P.Val']))
            logFC = np.array(list(result_py['logFC']))

        return self.sort(genes, p_values, logFC)