import pandas as pd
import logging
from pgmpy.estimators import PC
from networkx import DiGraph
import networkx as nx

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def generate_ccg_network(data, whitelist, blacklist, variable_to_tier, output_filename):
    """
        Constructs a Causal Conditional Graph (CCG) network based on the provided dataset. 
        It applies conditional independence tests to estimate the graph structure, 
        adjusts the graph by removing edges listed in a blacklist and adding edges from a whitelist, 
        and saves the network to a csv file.

        data: The dataset to be used to construct the CCG network.
        whitelist: List of edges (as tuples of node names) that should be added to the graph if not already present, based on prior knowledge or specific hypotheses.
        blacklist: List of edges (as tuples of node names) that should be removed from the graph to avoid spurious connections or known inaccuracies.
        variable_to_tier: A mapping from variable names to their respective tier levels, used to ensure the temporal direction of causality in the graph.
        output_filename: The path to the file where the network graph will be saved as a CSV.
    """
    network = DiGraph()
    try:
        logger.info("Processing combined dataset")
        pc = PC(data)
        pdag = pc.estimate(variant='stable', ci_test='chi_square', return_type='cpdag', significance_level=0.01, max_cond_vars=5)

        for edge in blacklist:
            if pdag.has_edge(*edge):
                pdag.remove_edge(*edge)

        for edge in whitelist:
            if not pdag.has_edge(*edge):
                pdag.add_edge(*edge)

        for u, v in pdag.edges():
            if variable_to_tier.get(u, 0) <= variable_to_tier.get(v, 0):
                network.add_edge(u, v)
            else:
                logger.info(f"Skipping edge from {u} to {v} as it would point backwards in time.")
                
        nx.write_edgelist(network, output_filename, data=False)
        logger.info(f"Network saved to {output_filename}")
    except Exception as e:
        logger.error(f"Error processing dataset: {e}")

def load_data(wave1_path, wave5_path):
    """
      Load and merge the datasets from two waves of a study. It reads data from CSV files for wave 1 and wave 5, 
      merges them based on the participant id, and cleans up ununsed column.
    """
    data_w1 = pd.read_csv(wave1_path)
    data_w5 = pd.read_csv(wave5_path)
    data = pd.merge(data_w1, data_w5, on='pident', how='inner')
    data = data.drop(columns=['pident', 'eincom01'])
    return data

if __name__ == "__main__":
    wave1_path = '../data/network/discrete_data_wave1.csv'
    wave5_path = '../data/network/discrete_data_wave5.csv'
    data = load_data(wave1_path, wave5_path)
    
    whitelist = [
        # Baseline known connections
        ('Sexe', 'aSerum_TG'), ('Sexe', 'aHDL_C'), ('Sexe', 'aApoB'),
        ('Sexe', 'asmokstat'), ('Sexe', 'acidep09'),
        ('asmokstat', 'aSerum_TG'), ('asmokstat', 'aHDL_C'), ('asmokstat', 'aApoB_'),
        ('aSerum_TG', 'amet_syn2'), ('aHDL_C', 'amet_syn2'), ('aApoB', 'amet_syn2'),
        ('ams_waist', 'amet_syn2'),
        ('ahsCRP', 'aIL6'),
        ('aSerum_TG', 'atri_med'), ('aHDL_C', 'ahdl_med'), ('aApoB', 'atri_med'),
        ('ams_gluc2', 'agluc_med'),
        ('ams_hpt', 'asbp_med'), ('ams_hpt', 'adbp_med'),
        ('acidep09', 'aedu'),
        ('aauditsc', 'asmokstat'),
        # Follow-up known connections
        ('aSerum_TG', 'eSerumTG'), ('aHDL_C', 'eHDLC'), ('aApoB', 'eApoB'),
        ('aIL6', 'eIL6'), ('ahsCRP', 'aHSCRP'),
        ('amet_syn2', 'emet_syn2'),
        ('ams_waist', 'ems_waist'),
        ('ams_hpt', 'ems_hpt'),
        ('ams_gluc2', 'ems_gluc2'),
        ('atri_med', 'etri_med'), ('ahdl_med', 'ehdl_med'), ('asbp_med', 'esbp_med'), ('adbp_med', 'edbp_med'), ('agluc_med', 'egluc_med'),
        ('acidep09', 'ecidep09'),
        ('aauditsc', 'eauditsc'),
        # Known connections between baseline and follow-up
        ('asmokstat', 'eSerum_TG'), ('asmokstat', 'eHDLC'), ('asmokstat', 'eApoB'),
        ('ahsCRP', 'eIL6'),
        ('acidep09', 'emet_syn2'),
        ('aedu', 'eauditsc')
    ]
    blacklist = [] 
    all_vars = list(data.columns)
    immutable_vars = ['Age', 'eage', 'Sexe', 'sex']
    for var in all_vars:
        for immutable_var in immutable_vars:
            if var != immutable_var:  # we don't blacklist a variable from causing itself
                blacklist.append((var, immutable_var))
    
    tiers = {
        'Context B': ['Age', 'Sexe'],
        'B': ['aedu', 'asmokstat', 'AIPMETO2', 'aauditsc', 'aIRSsum9', 'abaiscal', 'aids', 'acidep09', 'amet_syn2', 'ams_waist', 'ams_hpt', 'ams_trig2', 'ams_hdl2', 'ams_gluc2', 'atri_med', 'ahdl_med', 'asbp_med', 'adbp_med', 'agluc_med', 'ahsCRP', 'aIL6', 'aApoB', 'aHDL_C', 'aTotFA', 'aSerum_TG', 'aGp', 'aIle'],
        'Context FU': ['eage', 'sex'],
        'FU': ['eipmeto2', 'eauditsc', 'eIRSsum9', 'ebaiscal', 'eids', 'ecidep09', 'emet_syn2', 'ems_waist', 'ems_hpt', 'ems_trig2', 'ems_hdl2', 'ems_gluc2',
        'etri_med', 'ehdl_med', 'esbp_med', 'edbp_med', 'egluc_med', 'eApoB', 'eHDLC', 'eTotFA', 'eSerumTG', 'eGp', 'eIle', 'eHSCRP', 'eIL6',]
    }
    
    tier_levels = {'Context B': 1, 'B': 2, 'Context FU': 3, 'FU': 4}

    # Mapping variables to their respective tier levels
    variable_to_tier = {}
    for tier_name, variables in tiers.items():
        tier_level = tier_levels[tier_name]
        for var in variables:
            variable_to_tier[var] = tier_level

    generate_ccg_network(data, whitelist, blacklist, variable_to_tier, "network.csv")
