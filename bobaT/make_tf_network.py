"""
Network construction and feature selection utilities for bobaT.net module.
 
This module provides functions for:
- Building gene regulatory networks from DIRECT-NET and FIGR data
- LASSO-based feature selection for network pruning
- Network comparison and analysis
"""
from . import enrichr
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split, GridSearchCV, KFold
import time

def prune(G_orig, prune_sources=True, prune_sinks=True):
    """Prune graph to remove nodes with no incoming or outgoing edges

    :param G_orig: NetworkX graph to prune
    :type G_orig: networkx.DiGraph
    :param prune_sources: Remove nodes with no incoming edges, defaults to True
    :type prune_sources: bool, optional
    :param prune_sinks: Remove nodes with no outgoing edges, defaults to True
    :type prune_sinks: bool, optional
    :return: Pruned network
    :rtype: networkx.DiGraph
    """
    G = G_orig.copy()
    n = len(G.nodes())
    nold = n + 1

    while n != nold:
        nold = n
        for tf in list(G.nodes()):
            if prune_sources == True:
                if G.in_degree(tf) == 0:
                    G.remove_node(tf)
            if prune_sinks == True:
                if G.out_degree(tf) == 0:
                    G.remove_node(tf)
            else:
                if G.in_degree(tf) == 0 and G.out_degree(tf) == 0:
                    G.remove_node(tf)
        n = len(G.nodes())
    return G


def prune_info(G_orig, prune_self_loops=True):
    """Prune graph to only include edges with evidence in multiple databases

    :param G_orig: NetworkX graph to prune
    :type G_orig: networkx.DiGraph
    :param prune_self_loops: Whether to prune self-loops in the network (edges where parent node == child node), defaults to True
    :type prune_self_loops: bool, optional
    :return: Pruned network
    :rtype: networkx.DiGraph
    """
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if "db" not in edges[target]:
                G.remove_edge(tf, target)
            elif len(edges[target]["db"]) < 2:
                G.remove_edge(tf, target)
    return prune(G)


def prune_to_chea(G_orig, prune_self_loops=True):
    """Prune graph to only include edges with evidence in ChEA databases

    :param G_orig: NetworkX graph to prune
    :type G_orig: networkx.DiGraph
    :param prune_self_loops: Whether to prune self-loops in the network (edges where parent node == child node), defaults to True
    :type prune_self_loops: bool, optional
    :return: Pruned network
    :rtype: networkx.DiGraph
    """
    G = G_orig.copy()
    for tf in list(G.nodes()):
        edges = G.adj[tf]
        for target in list(edges.keys()):
            if tf == target and prune_self_loops:
                G.remove_edge(tf, target)
                continue
            if "db" in edges[target]:
                if not True in ["ChEA" in i for i in edges[target]["db"]]:
                    G.remove_edge(tf, target)
    #                if len(edges[target]['db']) < 2: G.remove_edge(tf, target)
    return prune(G)


def make_network(
    tfs,
    outdir="",
    do_prune=True,
    prune_sinks=True,
    prune_sources=True,
    do_prune_info=True,
    prune_self_loops=True,
    do_prune_to_chea=True,
    save_unfiltered=False,
    network_name="network",
):
    """Make network from list of tfs and save as csv files with various levels of pruning.

    :param tfs: List of transcription factor gene names that will be searched in enrichR databases.
    :type tfs: List[str]
    :param outdir: Output directory where network csvs will be saved, defaults to ""
    :type outdir: str or None
    :param do_prune: Prune network to remove nodes with no sinks and/or sources and generate a new network file called <network_name>_pruned.csv, defaults to True
    :type do_prune: bool, optional
    :param prune_sinks: Whether to prune sink nodes from the network (no outgoing nodes), defaults to True
    :type prune_sinks: bool, optional
    :param prune_sources: Whether to prune source nodes the network (no incoming nodes), defaults to True
    :type prune_sources: bool, optional
    :param do_prune_info: Whether to prune the network to edges with evidence in more than one database from enrichR and generate a new network file called <network_name>_high_evidence.csv, defaults to True
    :type do_prune_info: bool, optional
    :param prune_self_loops: Whether to prune self-loops in the network (edges where parent node == child node), defaults to True
    :type prune_self_loops: bool, optional
    :param do_prune_to_chea: Whether to prune the network to only edges with evidence in ChEA databases and generate a new network file called <network_name>_chea.csv, defaults to True
    :type do_prune_to_chea: bool, optional
    :param save_unfiltered: Whether to save the unfiltered network before pruning, defaults to False. Note that if all pruning options are False and this is set to False, no network will be saved.
    :type save_unfiltered: bool, optional
    :param network_name: Prefix name of network csv files, defaults to `network`
    :type network_name: str
    :return: Network with highest level of pruning
    :rtype: NetworkX Graph
    """

    G = nx.DiGraph()
    # prelim_G = nx.DiGraph()
    # with open("/Users/sarahmaddox/Dropbox (Vanderbilt)/Quaranta_Lab/SCLC/Network/mothers_network.csv") as infile:
    #     for line in infile:
    #         line = line.strip().split(',')
    #         prelim_G.add_edge(line[0], line[1])

    for tf in tfs:
        G.add_node(tf)

    for tf in tfs:
        enrichr.build_tf_network(G, tf, tfs)
        time.sleep(1)

    # for edge in prelim_G.edges():
    #     if edge[0] in tfs and edge[1] in tfs:
    #         G.add_edge(edge[0], edge[1])

    if save_unfiltered:
        outfile = open(
            os.path.join(outdir, f"{network_name}_unfiltered.csv"),
            "w",
        )
        for edge in G.edges():
            outfile.write("%s,%s\n" % (edge[0], edge[1]))
        outfile.close()

    if do_prune:
        Gp = prune(G, prune_sinks=prune_sinks, prune_sources=prune_sources)

        outfile = open(
            os.path.join(outdir, f"{network_name}_pruned.csv"),
            "w",
        )
        for edge in Gp.edges():
            outfile.write("%s,%s\n" % (edge[0], edge[1]))
        outfile.close()
    else:
        Gp = G

    if do_prune_info:
        Gpp = prune_info(Gp, prune_self_loops=prune_self_loops)

        outfile = open(
            os.path.join(outdir, f"{network_name}_high_evidence.csv"),
            "w",
        )
        for edge in Gpp.edges():
            outfile.write("%s,%s\n" % (edge[0], edge[1]))
        outfile.close()
    else:
        Gpp = Gp

    if do_prune_to_chea:
        Gppp = prune_to_chea(Gpp, prune_self_loops=prune_self_loops)

        outfile = open(
            os.path.join(outdir, f"{network_name}_chea.csv"),
            "w",
        )
        for edge in Gppp.edges():
            outfile.write("%s,%s\n" % (edge[0], edge[1]))
        outfile.close()
    else:
        Gppp = Gpp

    return Gppp


 
 
def save_network(network_file, G, attributes=True, overwrite=True):
    """
    Save a NetworkX graph to a CSV file.
    
    Parameters
    ----------
    network_file : str
        Path to output file
    G : networkx.DiGraph
        Network graph to save
    attributes : bool, default=True
        Whether to save edge attributes
    overwrite : bool, default=True
        Whether to overwrite existing file
        
    Returns
    -------
    bool
        True if saved, False if file exists and overwrite=False
    """
    from os.path import exists
    
    if exists(network_file) and not overwrite:
        print(f"Network file {network_file} already exists. Set overwrite=True to replace.")
        return False
    
    with open(network_file, "w") as outfile:
        for edge in G.edges():
            outfile.write(f"{edge[0]},{edge[1]}")
            if attributes:
                for attr in G[edge[0]][edge[1]].keys():
                    outfile.write(f",{G[edge[0]][edge[1]][attr]}")
            outfile.write("\n")
    
    print(f"Network saved to {network_file}")
    return True
 
 
def add_connections(G, net_df, threshold_var=None, threshold=0, 
                   parent_node='TF motif', child_node='Target_gene',
                   add_weight=True, weight='motif score', evidence='direct-net'):
    """
    Add edges to a NetworkX graph from a dataframe.
    
    Parameters
    ----------
    G : networkx.DiGraph
        Network graph to add edges to
    net_df : pd.DataFrame
        DataFrame containing edge information
    threshold_var : str, optional
        Column name to threshold on
    threshold : float, default=0
        Threshold value for filtering edges
    parent_node : str, default='TF motif'
        Column name for source nodes
    child_node : str, default='Target_gene'
        Column name for target nodes
    add_weight : bool, default=True
        Whether to add edge weights
    weight : str, default='motif score'
        Column name for edge weights
    evidence : str, default='direct-net'
        Evidence type for edge annotation
        
    Returns
    -------
    networkx.DiGraph
        Updated graph with new edges
    """
    for i, r in net_df.iterrows():
        if G.has_edge(r[parent_node], r[child_node]):
            continue
        
        # Check threshold if specified
        if threshold_var is not None:
            if r[threshold_var] <= threshold:
                continue
        
        # Add edge with or without weight
        if add_weight:
            G.add_edge(r[parent_node], r[child_node], 
                      weight=r[weight], evidence=evidence)
        else:
            G.add_edge(r[parent_node], r[child_node], evidence=evidence)
    
    return G
 
 
def build_network_from_directnet(direct_net_df, threshold=0, threshold_var='motif score'):
    """
    Build a gene regulatory network from DIRECT-NET data.
    
    Parameters
    ----------
    direct_net_df : pd.DataFrame
        DIRECT-NET data with columns: 'TF motif', 'Target_gene', 'motif score'
    threshold : float, default=0
        Minimum motif score threshold
    threshold_var : str, default='motif score'
        Column name to threshold on
        
    Returns
    -------
    networkx.DiGraph
        Gene regulatory network
    """
    # Get all unique TFs and targets
    tfs = list(set(direct_net_df['TF motif'].tolist() + 
                   direct_net_df['Target_gene'].tolist()))
    
    # Initialize directed graph
    G = nx.DiGraph()
    G.add_nodes_from(tfs)
    
    # Add connections
    G = add_connections(G, direct_net_df, threshold=threshold, 
                       threshold_var=threshold_var)
    
    print(f"Network built with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G
 
 
def lasso_feature_selection(network, data, network_name, output_dir, 
                            alphas=None, save_network=True, plot=True):
    """
    Perform LASSO regression for feature selection on network edges.
    
    For each target gene, use LASSO to identify the most important 
    transcription factor regulators from the network.
    
    Parameters
    ----------
    network : pd.DataFrame
        Network dataframe with columns: 'source', 'target', 'weight', 'evidence'
    data : pd.DataFrame
        Expression data with genes as columns
    network_name : str
        Name for output files
    output_dir : str
        Directory for saving results
    alphas : list, optional
        Alpha values to test. Default: [0.00001, 0.0001, 0.001, 0.01, 0.1, 1]
    save_network : bool, default=True
        Whether to save filtered networks
    plot : bool, default=True
        Whether to generate coefficient plots
        
    Returns
    -------
    dict
        Dictionary mapping alpha values to best scores
    """
    if alphas is None:
        alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1]
    
    # Create output directory
    base_dir = os.path.join(output_dir, "feature_selection", network_name)
    os.makedirs(base_dir, exist_ok=True)
    
    targets = list(np.unique(network['target']))
    alpha_scores = {alpha: [] for alpha in alphas}
    
    print(f"Running LASSO feature selection for {len(targets)} targets...")
    
    for target in targets:
        # Get TFs for this target
        self_target = False
        tfs = list(np.unique(network[network['target'] == target]['source']))
        
        if target in tfs:
            self_target = True
            tfs.remove(target)
        
        if len(tfs) == 0:
            print(f"Skipping {target}: no regulators found")
            continue
        
        # Prepare data
        X = data[tfs]
        y = data[target]
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=0
        )
        
        # Grid search for best alpha
        params = {"alpha": np.linspace(0.00001, 1, 50)}
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        lasso = Lasso()
        lasso_cv = GridSearchCV(lasso, param_grid=params, cv=kf)
        lasso_cv.fit(X, y)
        
        best_alpha = lasso_cv.best_params_['alpha']
        best_score = lasso_cv.best_score_
        
        # Test each specified alpha
        for alpha in alphas:
            feature_sel_folder = os.path.join(base_dir, f"Lasso_alpha_{alpha}")
            os.makedirs(feature_sel_folder, exist_ok=True)
            
            lasso_model = Lasso(alpha=alpha)
            lasso_model.fit(X_train, y_train)
            score = lasso_model.score(X_test, y_test)
            alpha_scores[alpha].append(score)
            
            # Get absolute coefficients
            lasso_coef = np.abs(lasso_model.coef_)
            
            if plot:
                plt.figure(figsize=(10, 6))
                plt.bar(tfs, lasso_coef)
                plt.xticks(rotation=90)
                plt.grid(alpha=0.3)
                plt.title(
                    f"Feature Selection: {target} (α={alpha})\n"
                    f"Score: {score:.3f} "
                    f"({score/best_score*100:.1f}% of best)"
                )
                plt.xlabel("Transcription Factors")
                plt.ylabel("Absolute Coefficient")
                plt.tight_layout()
                
                filename = f"{target}_best.png" if alpha == best_alpha else f"{target}.png"
                plt.savefig(os.path.join(feature_sel_folder, filename))
                plt.close()
            
            # Save filtered network
            if save_network:
                feature_subset = np.array(tfs)[lasso_coef > 0]
                if self_target:
                    feature_subset = np.append(feature_subset, target)
                
                subset = network[
                    (network['target'] == target) & 
                    (network['source'].isin(feature_subset))
                ]
                
                output_file = os.path.join(
                    feature_sel_folder, 
                    f"{network_name}_Lasso_{alpha}.csv"
                )
                
                # Append to file or create new
                if os.path.isfile(output_file):
                    subset.to_csv(output_file, mode='a', header=False, index=False)
                else:
                    subset.to_csv(output_file, index=False)
    
    print("LASSO feature selection complete!")
    return alpha_scores
 
 
def compare_networks(original_network, filtered_network, alpha, 
                    save=True, save_dir=""):
    """
    Compare original network to filtered network and visualize differences.
    
    Parameters
    ----------
    original_network : pd.DataFrame
        Original network with 'source' and 'target' columns
    filtered_network : pd.DataFrame
        Filtered network with 'source' and 'target' columns
    alpha : float
        Alpha value used for filtering
    save : bool, default=True
        Whether to save the histogram
    save_dir : str, default=""
        Directory for saving plots
        
    Returns
    -------
    dict
        Summary statistics about the network comparison
    """
    # Merge to find differences
    both = original_network.merge(
        filtered_network, 
        on=['source', 'target'],
        how='left', 
        indicator=True
    )
    difference = both[both["_merge"] != "both"]
    
    # Calculate statistics
    n_dropped = difference.shape[0]
    targets = list(np.unique(filtered_network['target']))
    
    parent_counts = []
    max_parents = 0
    max_parents_target = ""
    max_parents_tfs = []
    
    for target in targets:
        tmp = filtered_network[filtered_network['target'] == target]
        n_parents = len(np.unique(tmp['source']))
        parent_counts.append(n_parents)
        
        if n_parents > max_parents:
            max_parents = n_parents
            max_parents_target = target
            max_parents_tfs = list(np.unique(tmp['source']))
    
    # Print summary
    print(f"\nNetwork comparison (α={alpha}):")
    print(f"  Dropped edges: {n_dropped}")
    print(f"  Target with most parents: {max_parents_target} ({max_parents} parents)")
    print(f"    Parents: {max_parents_tfs}")
    
    # Plot histogram
    bins = np.arange(max_parents + 2) - 0.5
    plt.figure(figsize=(10, 6))
    plt.hist(parent_counts, bins=bins, edgecolor='#e0e0e0')
    plt.xticks(range(max_parents + 1))
    plt.xlim([-1, max_parents + 1])
    plt.xlabel("Number of Regulators")
    plt.ylabel("Number of Target Genes")
    plt.title(f"Regulator Distribution (α={alpha})")
    plt.tight_layout()
    
    if save and save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, f"parent_hist_{alpha}.png"))
        plt.close()
    else:
        plt.show()
    
    return {
        'alpha': alpha,
        'n_dropped': n_dropped,
        'max_parents': max_parents,
        'max_parents_target': max_parents_target,
        'parent_counts': parent_counts
    }
 
 
def adaptive_lasso_pruning(network, output_dir, network_name, 
                          max_parents=8, alphas=None):
    """
    Adaptively select alpha for each target to achieve desired max parents.
    
    For targets with more than max_parents regulators, find the smallest
    alpha that reduces regulators to max_parents or fewer.
    
    Parameters
    ----------
    network : pd.DataFrame
        Network dataframe with 'source' and 'target' columns
    output_dir : str
        Directory containing LASSO results
    network_name : str
        Name of the network
    max_parents : int, default=8
        Maximum number of regulators per target
    alphas : list, optional
        Alpha values to search. Default: [0.00001, 0.0001, 0.001, 0.01, 0.1, 1]
        
    Returns
    -------
    pd.DataFrame
        Pruned network with at most max_parents regulators per target
    dict
        Dictionary mapping each target to its selected alpha
    """
    if alphas is None:
        alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1]
    
    tfs = sorted(list(set(network['source']).union(set(network['target']))))
    tf_dict = {}
    tf_parent_dict = {}
    parent_counts = []
    
    base_dir = os.path.join(output_dir, "feature_selection", network_name)
    
    print(f"Adaptively pruning network to max {max_parents} parents per target...")
    
    for tf in tfs:
        tmp = network[network['target'] == tf]
        n_parents = len(np.unique(tmp['source']))
        
        # If already below threshold, keep all parents
        if n_parents <= max_parents:
            tf_dict[tf] = 0
            tf_parent_dict[tf] = list(np.unique(tmp['source']))
            parent_counts.append(n_parents)
            continue
        
        # Search for appropriate alpha
        found = False
        for alpha in alphas:
            network_file = os.path.join(
                base_dir, 
                f"Lasso_alpha_{alpha}",
                f"{network_name}_Lasso_{alpha}.csv"
            )
            
            if not os.path.exists(network_file):
                print(f"Warning: {network_file} not found")
                continue
            
            filtered_net = pd.read_csv(network_file)
            tmp_filtered = filtered_net[filtered_net['target'] == tf]
            n_parents_filtered = len(np.unique(tmp_filtered['source']))
            
            if n_parents_filtered <= max_parents:
                tf_dict[tf] = alpha
                tf_parent_dict[tf] = list(np.unique(tmp_filtered['source']))
                parent_counts.append(n_parents_filtered)
                found = True
                break
        
        if not found:
            print(f"Warning: Could not reduce {tf} to {max_parents} parents")
            tf_dict[tf] = "NA"
            tf_parent_dict[tf] = []
            parent_counts.append(0)
    
    # Create pruned network dataframe
    pruned_edges = []
    for target, sources in tf_parent_dict.items():
        for source in sources:
            pruned_edges.append({'source': source, 'target': target})
    
    pruned_network = pd.DataFrame(pruned_edges)
    
    # Plot distribution
    max_parents_final = max(parent_counts) if parent_counts else 0
    bins = np.arange(max_parents_final + 2) - 0.5
    
    plt.figure(figsize=(10, 6))
    plt.hist(parent_counts, bins=bins, edgecolor='#e0e0e0')
    plt.xticks(range(max_parents_final + 1))
    plt.xlim([-1, max_parents_final + 1])
    plt.xlabel("Number of Regulators")
    plt.ylabel("Number of Target Genes")
    plt.title(f"Regulator Distribution (Adaptive LASSO, max={max_parents})")
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, f"parent_hist_adaptive_max{max_parents}.png"))
    plt.close()
    
    print(f"Adaptive pruning complete! Final network: {len(pruned_network)} edges")
    
    return pruned_network, tf_dict
 