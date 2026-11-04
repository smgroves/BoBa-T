import json
import requests
import numpy as np
import networkx as nx
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
import pandas as pd
from matplotlib import pyplot as plt
import time
from requests.exceptions import ConnectionError

def get_libraries():
    return """
[u'Kinase_Perturbations_from_L1000',
 u'Jensen_COMPARTMENTS',
 u'ENCODE_TF_ChIP-seq_2015',
 u'ENCODE_TF_ChIP-seq_2014',
 u'TargetScan_microRNA',
 u'ESCAPE',
 u'MCF7_Perturbations_from_GEO_up',
 u'LINCS_L1000_Kinase_Perturbations_down',
 u'Jensen_DISEASES',
 u'MGI_Mammalian_Phenotype_Level_4',
 u'MGI_Mammalian_Phenotype_Level_3',
 u'Drug_Perturbations_from_GEO_up',
 u'NCI-Nature_2016',
 u'NCI-Nature_2015',
 u'GO_Molecular_Function_2013',
 u'Disease_Signatures_from_GEO_up_2014',
 u'Allen_Brain_Atlas_down',
 u'MSigDB_Computational',
 u'ENCODE_Histone_Modifications_2015',
 u'MGI_Mammalian_Phenotype_2017',
 u'Chromosome_Location',
 u'LINCS_L1000_Chem_Pert_down',
 u'Disease_Perturbations_from_GEO_up',
 u'Jensen_TISSUES',
 u'MGI_Mammalian_Phenotype_2013',
 u'Virus_Perturbations_from_GEO_up',
 u'LINCS_L1000_Kinase_Perturbations_up',
 u'Transcription_Factor_PPIs',
 u'Single_Gene_Perturbations_from_GEO_up',
 u'ChEA_2013',
 u'KEA_2013',
 u'ChEA_2016',
 u'ChEA_2015',
 u'Drug_Perturbations_from_GEO_2014',
 u'Drug_Perturbations_from_GEO_down',
 u'Aging_Perturbations_from_GEO_down',
 u'TRANSFAC_and_JASPAR_PWMs',
 u'Kinase_Perturbations_from_GEO',
 u'GO_Biological_Process_2015',
 u'NCI-60_Cancer_Cell_Lines',
 u'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
 u'LINCS_L1000_Chem_Pert_up',
 u'WikiPathways_2015',
 u'WikiPathways_2016',
 u'MSigDB_Oncogenic_Signatures',
 u'Genes_Associated_with_NIH_Grants',
 u'Kinase_Perturbations_from_GEO_down',
 u'WikiPathways_2013',
 u'PPI_Hub_Proteins',
 u'Cancer_Cell_Line_Encyclopedia',
 u'Disease_Signatures_from_GEO_down_2014',
 u'Human_Gene_Atlas',
 u'Pfam_InterPro_Domains',
 u'TF-LOF_Expression_from_GEO',
 u'GeneSigDB',
 u'Disease_Perturbations_from_GEO_down',
 u'Mouse_Gene_Atlas',
 u'GO_Molecular_Function_2015',
 u'LINCS_L1000_Ligand_Perturbations_down',
 u'Epigenomics_Roadmap_HM_ChIP-seq',
 u'LINCS_L1000_Ligand_Perturbations_up',
 u'KEA_2015',
 u'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
 u'Ligand_Perturbations_from_GEO_up',
 u'Single_Gene_Perturbations_from_GEO_down',
 u'Genome_Browser_PWMs',
 u'Ligand_Perturbations_from_GEO_down',
 u'Kinase_Perturbations_from_GEO_up',
 u'GO_Cellular_Component_2015',
 u'Allen_Brain_Atlas_up',
 u'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO',
 u'Human_Phenotype_Ontology',
 u'Microbe_Perturbations_from_GEO_up',
 u'HomoloGene',
 u'MCF7_Perturbations_from_GEO_down',
 u'Old_CMAP_up',
 u'ENCODE_Histone_Modifications_2013']
"""

# absolute max
def amax(a,b):
    aa = a
    bb = b
    if a < 0: aa*=-1
    if b < 0: bb*=-1
    if aa < bb: return b
    return a
    
# absolute min
def amin(a,b):
    aa = a
    bb = b
    if a < 0: aa*=-1
    if b < 0: bb*=-1
    if aa < bb: return a
    return b

def query_gene(gene):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/genemap'
    query_string = '?json=true&setup=true&gene=%s'%gene
    try:
        response = requests.get(ENRICHR_URL + query_string)
        if not response.ok:
            raise Exception('Error searching for terms')
        return json.loads(response.text)
    except ConnectionError as e:
        print("Connection rejected on %s, waiting 6 seconds and trying again"%gene)
        time.sleep(6)
        return query_gene(gene)

def submit_gene_list(genes, description = "Quaranta gene set TF analysis"):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join(genes)

    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
#    return data
    return data['userListId']

def enrich(user_list_id, library):
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    
    response = requests.get(ENRICHR_URL + query_string % (user_list_id, library))

    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)[library]
    # [ index, target (somewhere in the string), p-value, z-score, combined score, [TFs], adjusted p-value]
    return data

def process_background(data, background, G = None):
    if G is None: G = nx.DiGraph()
    for d in data:
        if ("ChEA" in background): source = get_chea_source(d)
        elif ("TRANSFAC" in background): source = get_transfac_source(d)
        else: raise Exception("Unemplemented background")
        targets = [str(i) for i in d[5]]
        p_value = d[2]
        z_score = d[3]
        combined_score = d[4]
        adj_p = d[6]
        if not G.has_node(source): G.add_node(source)
        for target in targets:
            if G.adj[source].has_key(target):
                G.adj[source][target]['p_value'] = min(G.adj[source][target]['p_value'], p_value)
                G.adj[source][target]['z_score'] = amax(G.adj[source][target]['z_score'], z_score)
                G.adj[source][target]['combined_score'] = max(G.adj[source][target]['combined_score'], combined_score)
                G.adj[source][target]['adj_p'] = min(G.adj[source][target]['adj_p'], adj_p)
            else:
                G.add_edge(source,target,p_value = p_value, z_score = z_score, combined_score = combined_score, adj_p = adj_p)
    return G

def process_chea_lists(data, G = None):
    if G is None: G = nx.DiGraph()
    for d in data:
        source = str(d[1].split('_')[0])
        targets = [str(i) for i in d[5]]
        p_value = d[2]
        z_score = d[3]
        combined_score = d[4]
        adj_p = d[6]
        if not G.has_node(source): G.add_node(source)
        for target in targets:
            if G.adj[source].has_key(target):
                G.adj[source][target]['p_value'].append(p_value)
                G.adj[source][target]['z_score'].append(z_score)
                G.adj[source][target]['combined_score'].append(combined_score)
                G.adj[source][target]['adj_p'].append(adj_p)
            else:
                G.add_edge(source,target,p_value = [p_value,], z_score = [z_score,], combined_score = [combined_score,], adj_p = [adj_p,])
    return G


# TODO: Deal with aliases!!!!!!!!!!!!!
def _first_token(term):
    return [term.split(' ')[0]]

def _jaspar_2025_sources(term):
    # dimer motifs are named "TF1::TF2"; credit each partner
    return term.split(' ')[0].split('::')

# Parsers for the TF name(s) in each library's term names, matching current Enrichr formats.
# In each of these libraries a term is a TF and its gene set is that TF's targets, so a hit
# for `gene` in term `t` means an edge source(t) -> gene.
LIBRARY_PARSERS = {
    'ChEA_2013': lambda t: [t.split('-')[0]],                      # 'EZH2-18974828-MESC-MOUSE'
    'ChEA_2015': _first_token,                                     # 'SUZ12 18555785 ChIP-Seq MESC Mouse'
    'ChEA_2016': _first_token,
    'ChEA_2022': _first_token,
    'ENCODE_TF_ChIP-seq_2014': _first_token,                       # 'CTCF H1'
    'ENCODE_TF_ChIP-seq_2015': _first_token,                       # 'MAX H1-hESC hg19'
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X': _first_token,     # 'EGR1 CHEA'
    'TRANSFAC_and_JASPAR_PWMs': _first_token,                      # 'TEAD2 (human)'
    'JASPAR_PWM_Human_2025': _jaspar_2025_sources,                 # 'ELK1::HOXB13'
    'JASPAR_PWM_Mouse_2025': _jaspar_2025_sources,                 # 'NEUROD2'
    'TRRUST_Transcription_Factors_2019': _first_token,             # 'FOXA2 mouse'
    'TF_Perturbations_Followed_by_Expression': _first_token,       # 'SOX3 KO MOUSE GSE53760 ...'
    'Rummagene_transcription_factors': lambda t: [t.split(' ')[-1]],  # '...-Gene Name AR'
}

# Binding (ChIP), motif, and curated TF-target libraries available in Enrichr as of 2026.
# TF_Perturbations_Followed_by_Expression (indirect, expression-based) and
# Rummagene_transcription_factors (mined supplementary tables) are parsed but not included;
# pass them explicitly in `libraries` to use them.
CURRENT_TF_LIBRARIES = [
    'ChEA_2013', 'ChEA_2015', 'ChEA_2016', 'ChEA_2022',
    'ENCODE_TF_ChIP-seq_2014', 'ENCODE_TF_ChIP-seq_2015',
    'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
    'TRANSFAC_and_JASPAR_PWMs', 'JASPAR_PWM_Human_2025', 'JASPAR_PWM_Mouse_2025',
    'TRRUST_Transcription_Factors_2019',
]

def build_tf_network(G, tf, tfs, libraries='current'):
    """Add edges into `tf` from any TF in `tfs` that Enrichr lists as a regulator of it.

    :param G: graph to add edges to; each edge stores its supporting libraries in ``db``
    :type G: networkx.DiGraph
    :param tf: gene to query
    :type tf: str
    :param tfs: candidate regulators; only edges whose source is in this list are added
    :type tfs: List[str]
    :param libraries: ``"current"`` (default) uses ``CURRENT_TF_LIBRARIES``; a list of
        names from ``LIBRARY_PARSERS`` uses those. Both use parsers matching current
        Enrichr term names. ``"legacy"`` (or None) keeps the original library list and
        term parsing; with current Enrichr term names, only ChEA_2013 and
        TRANSFAC_and_JASPAR_PWMs yield edges in that mode.
    :type libraries: None, str, or List[str]
    :return: G with edges added
    :rtype: networkx.DiGraph
    """
    if libraries is not None and libraries != 'legacy':
        if libraries == 'current': libraries = CURRENT_TF_LIBRARIES
        unknown = [lib for lib in libraries if lib not in LIBRARY_PARSERS]
        if unknown:
            raise ValueError("No parser for libraries %s; known: %s" % (unknown, sorted(LIBRARY_PARSERS)))
        print(tf)
        enrichr = query_gene(tf)
        for lib in libraries:
            if lib not in enrichr['gene']: continue
            sources = set()
            for term in enrichr['gene'][lib]:
                sources.update(str(s).upper() for s in LIBRARY_PARSERS[lib](term))
            for source in sources:
                if source in tfs:
                    if tf in G.adj[source]:
                        if lib not in G.adj[source][tf]['db']: G.adj[source][tf]['db'].append(lib)
                    else: G.add_edge(source, tf, db=[lib,])
        return G

#    G = nx.DiGraph()
#    for gene in tfs:
#        G.add_node(gene)
#    for gene in tfs:
    for gene in [tf,]:
#        time.sleep(3)
        print(gene)
        enrichr = query_gene(gene)
        if ('ChEA_2013' in enrichr['gene']):
            for source in {str(i.split('-')[0]) for i in enrichr['gene']['ChEA_2013']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append("ChEA_2013")
                    else: G.add_edge(source,gene,db=["ChEA_2013",])

        if ('ChEA_2015' in enrichr['gene']):
            for source in {str(i.split('_')[0]) for i in enrichr['gene']['ChEA_2015']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append("ChEA_2015")
                    else: G.add_edge(source,gene,db=["ChEA_2015",])
                
        if ('ChEA_2016' in enrichr['gene']):
            for source in {str(i.split('_')[0]) for i in enrichr['gene']['ChEA_2016']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append("ChEA_2016")
                    else: G.add_edge(source,gene,db=["ChEA_2016",])

        if ('ENCODE_TF_ChIP-seq_2014' in enrichr['gene']):
            for source in {str(i.split('_')[0]) for i in enrichr['gene']['ENCODE_TF_ChIP-seq_2014']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append("ENCODE_TF_ChIP-seq_2014")
                    else: G.add_edge(source,gene,db=["ENCODE_TF_ChIP-seq_2014",])
            
        if ('ENCODE_TF_ChIP-seq_2015' in enrichr['gene']):
            for source in {str(i.split('_')[0]) for i in enrichr['gene']['ENCODE_TF_ChIP-seq_2015']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append("ENCODE_TF_ChIP-seq_2015")
                    else: G.add_edge(source,gene,db=["ENCODE_TF_ChIP-seq_2015",])

        if ('TRANSFAC_and_JASPAR_PWMs' in enrichr['gene']):
            for source in {str(i.split(' ')[0]) for i in enrichr['gene']['TRANSFAC_and_JASPAR_PWMs']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append('TRANSFAC_and_JASPAR_PWMs')
                    else: G.add_edge(source,gene,db=['TRANSFAC_and_JASPAR_PWMs',])
                
        if ('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X' in enrichr['gene']):
            for source in {str(i.split('_')[0]) for i in enrichr['gene']['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X']}:
                if source in tfs:
                    if gene in G.adj[source]: G.adj[source][gene]['db'].append('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X')
                    else: G.add_edge(source,gene,db=['ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',])
    return G
        

def prune_weak_edges(G):
    scores = []
    for e in G.adjs():
        scores.append(G.adj[e[0]][e[1]]['combined_score'])
    x = np.zeros((len(scores),1))
    x[:,0] = np.asarray(scores)
    km = KMeans(n_clusters=2)
    clusters = km.fit_predict(x)

    km_clusts = {0:[],1:[]}
    for i,s in zip(clusters,scores):
        km_clusts[i].append(s)
    
    density = gaussian_kde(scores)
    xs = np.linspace(min(scores),max(scores),200)
    
    Gp = nx.DiGraph()
    
    c0_mean = np.mean(km_clusts[0])
    c1_mean = np.mean(km_clusts[1])
    if c0_mean < c1_mean: cindex = 1
    else: cindex = 0
    
    cutoff = (max(km_clusts[1-cindex]) + min(km_clusts[cindex])) / 2.
    
    for c,e in zip(clusters,G.adjs()):
        if c == cindex: Gp.add_edge(e[0],e[1])
    
    plt.plot(xs,density(xs))
    plt.plot([cutoff,cutoff],[0,max(density(xs))])
    plt.show()
    
    return Gp
    
    
def get_chea_source(d):
    return str(d[1].split('_')[0])

def get_transfac_source(d):
    return str(d[1].split(' ')[0])
    
def get_targetscan_source(d):
    if d is not None: raise Exception("TargetScan code not actually implemented yet... it can return multiple miRNAs on a single line")
    return str(d[1].split(',')[0])
    
def get_browser_source(d):
    d = d[1]
    if ("UNKNOWN" in d): return "?"
    gene = d.split('$')[1].split('_')[0]
    if ("ALPHA" in gene[-5:]): gene = gene.replace("ALPHA","A")
    if ("BETA" in gene[-4:]): gene = gene.replace("BETA","B")
    if ("GAMMA" in gene[-5:]): gene = gene.replace("GAMMA","G")
    if ("DELTA" in gene[-5:]): gene = gene.replace("DELTA","D")
    if ("KAPPA" in gene[-5:]): gene = gene.replace("KAPPA","K")
    return str(gene)
    
def get_encode_source(d):
    return str(d[1].split('_')[0])
