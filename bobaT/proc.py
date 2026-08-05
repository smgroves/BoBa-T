from . import utils as ut

import os
import numpy as np


def binarize_data(
    data,
    phenotype_labels=None,
    threshold=0.5,
    save=False,
    save_dir=None,
    fname="binarized_data",
):
    if phenotype_labels is None:
        binaries = set()
    else:
        binaries = dict()
        for c in phenotype_labels["class"].unique():
            binaries[c] = set()

    f = np.vectorize(lambda x: "0" if x < threshold else "1")

    for sample in data.index:
        b = ut.state2idx("".join(f(data.loc[sample])))

        if phenotype_labels is None:
            binaries.add(b)
        else:
            binaries[phenotype_labels.loc[sample, "class"]].add(b)

    if save == True:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        with open(save_dir + os.sep + fname + ".csv", "w+") as outfile:
            for k in binaries.keys():
                outfile.write(f"{k}: {binaries[k]}\n")

    return binaries

