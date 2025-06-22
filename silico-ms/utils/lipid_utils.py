import re 
from itertools import chain
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns 

"""
from .chemical_utils import ELEMENT_TO_MASS

LIPDS_TO_SMILES = {
    'PC': None,
    'LPC': None,
    'PE': None, 
    'LPE': None,
    'PI': None, 
    'LPI': None,
    'PG': None, 
    'LPG': None,
    'PS': None, 
    'LPS': None,
    'Cer': None,
    'SM': None,
}


class Lipid:
    def __init__(self, 
                 string: str):
        self.lipid_class, self.fas = self.shorthand_to_lipid(string)
        self.sn = False
    
    def shorthand_to_lipid(self, 
                           shorthand: str, 
                           style: str = 'pubchem'):
        results = re.split(r'[()]', shorthand)
        results = [re.split(r'[_/]', s) for s in results]
        results = list(chain(*results))
        results = [s for s in results if len(s) > 0 ]
        lipid_class = results[0]
        fa_strings = results[1:]
        fas = [self.shorthand_to_fa(fa_string) for fa_string in fa_strings]
        return lipid_class , fas
    
    def shorthand_to_fa(self, 
                        shorthand: str):
        key_list = ["carbon", "nums_db"]
        value_list = shorthand.split(r':')
        value_list = [int(value) for value in value_list]
        fa_dict = dict(zip(key_list, value_list))
        fa_dict["shorthand"] = shorthand
        fa_dict["sn"] = 0
        fa_dict["pos_db"] = np.zeros(50)
        fa_dict["formula"] = 'C' + str(value_list[0]) + \
                             'H' + str(2 * value_list[0] - 2 * value_list[1]) + \
                             'O' + str(2)
        return fa_dict
    
    def lipid_to_shorthand(self) -> str:
        fas = '_'.join([fa["shorthand"] for fa in self.fas])
        lipid = self.lipid_class + "({fas})".format(fas=fas)
        return lipid
"""


def pubchem_to_lipidmaps(s):
    s = re.sub(r'\(', ' ', s, count=1)
    s = re.sub(r'\)$', '', s, count=1)
    return s


def get_lipid_subclass(shorthand, style = "pubchem"):
    match style:
        case "pubchem":
            shorthand = pubchem_to_lipidmaps(shorthand)
        case _:
            pass
    lipid_subclass = re.split(' ', shorthand)[0]
    return lipid_subclass


def get_species_level_lipid(shorthand, style = "pubchem"):
    match style:
        case "pubchem":
            shorthand = pubchem_to_lipidmaps(shorthand)
        case _:
            pass  
    shorthand = re.sub(r'\(.*?\)', '', shorthand)
    lipid_species = re.split(' ', shorthand)[0]
    lipid_fatty_acid = re.split(' ', shorthand)[1]
    lipid_fatty_acid_list = re.split('_', lipid_fatty_acid)
    
    def get_num_carbon(fa):
        num_carbon = re.split(r'\:', fa)[0]
        num_carbon = re.sub(r'([a-zA-Z])*', '', num_carbon)
        return int(num_carbon)
    
    def get_num_double_bond(fa):
        num_double_bond = re.split(r'\:', fa)[1]
        return int(num_double_bond)
    
    num_carbon = sum([get_num_carbon(fa) for fa in lipid_fatty_acid_list])
    num_double_bond = sum([get_num_double_bond(fa) for fa in lipid_fatty_acid_list])
    
    lipid = lipid_species + " " + str(num_carbon) + ":" + str(num_double_bond)
    return lipid    
    

def get_chain_level_lipid(shorthand, style = "pubchem"):
    match style:
        case "pubchem":
            shorthand = pubchem_to_lipidmaps(shorthand)
        case _:
            pass
    chain_level_lipid = re.sub(r'\(.*?\)', '', shorthand)
    return chain_level_lipid


def get_lipid_fatty_acid(shorthand, style = "pubchem"):
    match style:
        case "pubchem":
            shorthand = pubchem_to_lipidmaps(shorthand)
        case _:
            pass
    lipid_fatty_acid = re.split(' ', shorthand)[1]
    lipid_fatty_acid_list = re.split('_', lipid_fatty_acid)
    lipid_fatty_acid_list = [re.sub(r'\(.*?\)', '', fa) for fa in lipid_fatty_acid_list]
    return lipid_fatty_acid_list


def get_fatty_acid_double_bond_pos(fatty_acid_shorthand):
    if re.search(r"[mdt]", fatty_acid_shorthand):
        return ''
    fatty_acid_shorthand = re.sub(r'P-|O-', '', fatty_acid_shorthand)
    fatty_acid_shorthand = re.sub(r'\(\dOH\)', '', fatty_acid_shorthand)
    total_carbon = re.split(r':', fatty_acid_shorthand)[0]
    
    double_bond_pos = re.findall(r"\((.*?)\)", fatty_acid_shorthand)
    if len(double_bond_pos) == 0:
        return ''
    else:
        double_bond_pos = double_bond_pos[0]
    double_bond_pos_list = re.split(r'\,', double_bond_pos)
    double_bond_pos_omega_list = [int(total_carbon)- int(pos) for pos in double_bond_pos_list]
    double_bond_pos_omega = min(double_bond_pos_omega_list)
    double_bond_pos_omega = "n-" + str(double_bond_pos_omega)
    return double_bond_pos_omega


def get_lipid_double_bond_pos(shorthand, style="pubchem"):
    match style:
        case "pubchem":
            shorthand = pubchem_to_lipidmaps(shorthand)
        case _:
            pass
    lipid_fatty_acid = re.split(' ', shorthand)[1]
    lipid_fatty_acid_list = re.split('_', lipid_fatty_acid)
    lipid_double_bond_pos_list = [get_fatty_acid_double_bond_pos(fa) for fa in lipid_fatty_acid_list]
    return lipid_double_bond_pos_list


def get_max_intensity_lipid_adduct(data):
    #df = pd.read_csv(file)
    df = pd.melt(data, id_vars=["name", "adduct"], var_name="sample", value_name="intensity")
    df1 = df.groupby(["name", "adduct"]).sum(["intensity"]).reset_index()
    df2 = df1.groupby(["name"]).max(["intensity"]).reset_index()
    df3 = pd.merge(df1, df2, how="inner", on=["name", "intensity"])
    df_final = pd.merge(data, df3, how="inner", on=["name", "adduct"])
    df_final = df_final.drop(["intensity"], axis=1)
    return df_final


def get_lipid_info(data):
    df_fatty_acid = data.copy()
    df_double_bond = data.copy()
    
    df_fatty_acid["fatty_acid"] = df_fatty_acid["name"].apply(get_lipid_fatty_acid)
    df_fatty_acid = df_fatty_acid.explode("fatty_acid")
    df_double_bond["double_bond"] = df_double_bond["name"].apply(get_lipid_double_bond_pos)
    df_double_bond = df_double_bond.explode("double_bond")
    
    df_total = df_fatty_acid.copy()
    df_total["subclass"] = df_total["name"].apply(get_lipid_subclass)
    #df_total["fatty_acid"] = df_total["fatty_acid"].apply(lambda x: "FA " + x)
    df_total["fatty_acid"] = df_total["fatty_acid"].apply(lambda x: "C" + x)
    df_total["double_bond"] = df_double_bond["double_bond"]
    df_total["species_level"] = df_total["name"].apply(get_species_level_lipid)
    df_total["chain_level"] = df_total["name"].apply(get_chain_level_lipid)
    
    return df_total
    
def get_composition_heatmap(data):
    data = data[data["double_bond"] != ""]
    data = data.groupby(["sample", "fatty_acid", "double_bond"]).sum(["intensity"]).reset_index()
    
    data['total_intensity'] = data.groupby(['sample', 'fatty_acid'])['intensity'].transform('sum')
    data['percentage'] = (data['intensity'] / data['total_intensity']) * 100
    average_percentage = data.groupby(["fatty_acid", "double_bond"])['percentage'].mean().reset_index()
    df_plot = average_percentage.pivot(index="fatty_acid", columns="double_bond", values="percentage")
    
    return df_plot


def get_fc_barplot(data, df_group, group1 = "HFD+Fru", group2 = "HFD"):
    data = data[data["double_bond"] != ""]
    data = data.groupby(["sample", "fatty_acid", "double_bond"]).sum(["intensity"]).reset_index()
    
    data["full_fatty_acid"] = data.apply(lambda x: x["fatty_acid"] + "(" + x["double_bond"] + ")", axis=1)    
    df_fc = pd.merge(data, df_group, how="inner", on=["sample"])
    df_avg = df_fc.groupby(["group", "full_fatty_acid"])["intensity"].mean().reset_index()
    
    df_avg_wide = df_avg.pivot(index="full_fatty_acid", columns="group", values="intensity").reset_index()
    df_avg_wide["FC"] = df_avg_wide[group1] / df_avg_wide[group2]
    df_avg_wide = df_avg_wide[df_avg_wide["FC"] != 0]
    df_avg_wide["log2FC"] =  np.log2(df_avg_wide["FC"])
    
    return df_avg_wide
    

def get_relative_composition_barplot(data):
    data = data[data["double_bond"] != ""]
    data = data.groupby(["sample", "chain_level", "double_bond"]).sum(["intensity"]).reset_index()
    data['total_intensity'] = data.groupby(['sample', 'chain_level'])['intensity'].transform('sum')
    data['percentage'] = (data['intensity'] / data['total_intensity']) * 100
    
    df_plot = data.groupby(["chain_level", "double_bond"])["percentage"].mean().reset_index()
    
    return df_plot


def get_DEA_volcanoplot(data, 
                        df_group,
                        fc_threshold = 1.2,
                        pval_threshold = 0.05):
    import math
    import numpy as np
    import pandas as pd
    from scipy.stats import ttest_ind
    import statsmodels.stats.multitest as multitest
    
    log2_fc_thres = math.log2(fc_threshold)
    #log10_pval_thres = -math.log10(pval_threshold)

    group_type = list(set(df_group["group"].to_list()))

    df_data_wide = data.pivot(index="name", columns="sample", values="intensity")
    data_standardized = (df_data_wide - df_data_wide.mean()) / df_data_wide.std()
    data_standardized = data_standardized.T
    
    group1_idx = df_group[df_group["group"] == group_type[0]].index.to_list()
    group2_idx = df_group[df_group["group"] == group_type[1]].index.to_list()
    group1_data = data_standardized.iloc[group1_idx, :]
    group2_data = data_standardized.iloc[group2_idx, :]

# 存储统计测试结果
    de_results = pd.DataFrame(index=df_data_wide.index, columns=['p_value', 'log2_fold_change'])

    for gene in df_data_wide.index:
        _, p_value = ttest_ind(group1_data[gene], group2_data[gene], equal_var=False)
        log2_fold_change = np.log2(group1_data[gene].mean() / group2_data[gene].mean())
        de_results.loc[gene, 'p_value'] = p_value
        de_results.loc[gene, 'log2_fold_change'] = log2_fold_change


    def get_regulation(row):
        if (row["p_value"] > pval_threshold) or abs(row["log2_fold_change"]) < log2_fc_thres:
            return "No Significant"
        else:
            if row["log2_fold_change"] < -log2_fc_thres:
                return "Down"
            else:
                return "Up"

    de_results["significant"] = de_results.apply(get_regulation, axis=1)
    _, padj, _, _ = multitest.multipletests(de_results["p_value"], alpha=0.05, method='fdr_bh')
    de_results["fdr"] = padj
    de_results["log10_fdr"] = -np.log10(padj.tolist())
    
    return de_results
    

def get_double_bond_percentage(data):
    data = data[data["double_bond"] != ""]
    data = data.groupby(["sample", "chain_level", "double_bond"]).sum(["intensity"]).reset_index()
    
    data['total_intensity'] = data.groupby(["sample", "chain_level"])['intensity'].transform('sum')
    data['percentage'] = (data['intensity'] / data['total_intensity']) * 100
    #df_percentage = data.groupby(["sample","fatty_acid", "double_bond"])['percentage']
    
    return data


def get_chain_level_percentage(data, 
                               df_group,
                               group1 = "HFD+Fru", 
                               group2 = "HFD",
                               double_bond = "n-9",
                               pval_threshold = 0.05):
    import pandas as pd
    from scipy.stats import ttest_ind
    import statsmodels.stats.multitest as multitest
    
    
    data = data[data["double_bond"] == double_bond]
    cols = ["sample", "chain_level", "percentage"]
    data = data[cols]

    df_data_wide = data.pivot(index="chain_level", columns="sample", values="percentage")
    #data_standardized = (df_data_wide - df_data_wide.mean()) / df_data_wide.std()
    #data_standardized = data_standardized.T
    data_standardized = df_data_wide.T
    
    group1_idx = df_group[df_group["group"] == group1].index.to_list()
    group2_idx = df_group[df_group["group"] == group2].index.to_list()
    group1_data = data_standardized.iloc[group1_idx, :]
    group2_data = data_standardized.iloc[group2_idx, :]

    #data_original =  df_data_wide.T
    #group1_data_original = data_original.iloc[group1_idx, :]
    #group2_data_original = data_original.iloc[group2_idx, :]
    
    de_results = pd.DataFrame(index=df_data_wide.index, columns=['p_value', group1, group2])

    for lipid in df_data_wide.index:
        _, p_value = ttest_ind(group1_data[lipid], group2_data[lipid], equal_var=False)
        de_results.loc[lipid, 'p_value'] = p_value
        de_results.loc[lipid, group1] = group1_data[lipid].mean()
        de_results.loc[lipid, group2] = group2_data[lipid].mean()
        #de_results.loc[lipid, group1] = group1_data_original[lipid].mean()
        #de_results.loc[lipid, group2] = group2_data_original[lipid].mean()
        

    de_results["significant"] = de_results["p_value"] < pval_threshold
    _, padj, _, _ = multitest.multipletests(de_results["p_value"], alpha=0.05, method='fdr_bh')
    de_results["fdr"] = padj
    
    return de_results.reset_index()