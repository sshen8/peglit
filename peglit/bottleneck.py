"""
Cluster sequences and output top-scoring sequence per cluster.
"""
import warnings
import random
import numpy as np
from scipy.cluster.hierarchy import ClusterWarning
from sklearn.cluster import AgglomerativeClustering as HAC
from Levenshtein import distance as levenshtein_distance
from . import constants

def apply_bottleneck(heap_scores, heap, bottleneck=None, seed=None, verbose=False):
    """
    Cluster sequences and output top-scoring sequence per cluster.
    """
    # Process arguments with defaults as needed
    random.seed(seed)
    if bottleneck is None:
        bottleneck = constants.DEFAULT_BOTTLENECK
    # Pick best, randomly tiebreak if needed
    def _pick_best(scores, choices):
        idx_maxed = np.where(scores == np.max(scores))[0]
        idx_chosen = random.choice(idx_maxed)
        return choices[idx_chosen]
    # Can just pick best output
    if bottleneck == 1 and not verbose:
        return [_pick_best(heap_scores, heap)]
    # Calculate features for each linker sequence i.e. edit distance to all other linker sequences
    features = np.zeros((len(heap), len(heap)), dtype=int)
    for i, seq_x in enumerate(heap):
        for j, seq_y in enumerate(heap):
            features[i, j] = levenshtein_distance(seq_x, seq_y)
    # Can just pick best output
    if bottleneck == 1 and verbose:
        return [_pick_best(heap_scores, heap)], features
    # Cluster linker sequences
    warnings.filterwarnings( # ignore warning - we know it's a distance matrix
        action="ignore", category=ClusterWarning,
        message="scipy.cluster: The symmetric non-negative hollow observation matrix looks "
                "suspiciously like an uncondensed distance matrix")
    clusters = HAC(n_clusters=bottleneck, linkage="complete").fit_predict(features)
    # Output highest-scoring linker sequence from each cluster
    output = []
    heap = np.array(heap)
    heap_scores_mean = np.mean(heap_scores, axis=1)
    for cluster_num in range(bottleneck):
        idx_cluster = clusters == cluster_num
        heap_cluster = heap[idx_cluster]
        cluster_scores = heap_scores_mean[idx_cluster]
        output.append(_pick_best(cluster_scores, heap_cluster))
    if verbose:
        return output, features
    return output
