import csv
import os

class AfiMap:

    def _do_HG_U133(self, filename):
        with open(filename, newline='') as csvfile:
            map_file_reader = csv.reader(csvfile,
                                         delimiter=',',
                                         quotechar='"',
                                         quoting=csv.QUOTE_ALL,
                                         skipinitialspace=True)
            first = True
            for idx, row in enumerate(map_file_reader):
                
                if row[0] == '#':
                    continue

                if first:
                    header = row
                    first = False
                    continue

                gene_symbol = row[14].split('//')[0].strip()
                if gene_symbol == "---":
                    gene_symbol = row[9]
                self.map[row[0]] = {
                    'gene_title':row[13],
                    'gene_symbol': gene_symbol,
                    'entrez':row[18],
                }

    def _do_GPL96(self, filename):
        with open(filename, newline='') as csvfile:
            map_file_reader = csv.reader(csvfile,
                                         delimiter='\t',
                                         quotechar='"',
                                         quoting=csv.QUOTE_ALL,
                                         skipinitialspace=True)
            first = True
            for idx, row in enumerate(map_file_reader):

                if row[0][0] == '#':
                    continue

                if first:
                    first = False
                    continue

                if row[0] not in self.map:
                    self.map[row[0]] = {
                        'gene_title':row[8],
                        'gene_symbol': row[8],
                        'entrez':row[8],
                    }

    def __init__(self, folder_name):
        hg_filename = os.path.join(folder_name, "HG-U133_Plus_2.na36.annot.csv")
        gpl_filename = os.path.join(folder_name, "GPL96-57554.txt")
        self.map = {}
        self._do_HG_U133(hg_filename)
        self._do_GPL96(gpl_filename)

       
    def convert_afi_list_to_gene_symbols(self, gene_afi_list):
        new_list = []
        for gene_afi in gene_afi_list:
            new_list.append(self.map[gene_afi]['gene_symbol'])
        return new_list