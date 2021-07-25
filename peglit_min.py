from math import prod
import random
import heapq
import numpy as np
from scipy.special import expit as sigmoid
from sklearn.cluster import AgglomerativeClustering as HAC
from Levenshtein import distance as levenshtein_distance
import RNA # ViennaRNA

BASE_SYMBOLS = {
    "A": ("A",), "C": ("C",), "G": ("G",), "T": ("T",), "U": ("T",),
    "W": ("A", "T"), "S": ("C", "G"), "M": ("A", "C"),
    "K": ("G", "T"), "R": ("A", "G"), "Y": ("C", "T"),
    "B": ("C", "G", "T"), "D": ("A", "G", "T"), "H": ("A", "C", "T"), "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T")}

def apply_filters(seq_pre, seq_linker, seq_post, ac_thresh, u_thresh, n_thresh):
    """
    Returns False if any filter is failed i.e. AC content < ac_thresh OR consecutive Us > u_thresh
    OR consecutive Ns > n_thresh. Otherwise, True if all filters are passed. All thresholds have
    units nt (i.e. ac_thresh is not a percent). Ts are treated as Us.
    """
    # AC content
    if seq_linker.count("A") + seq_linker.count("C") < ac_thresh:
        return False
    # Consecutive U
    seq_neighborhood = seq_pre[-(u_thresh):] + seq_linker + seq_post[:u_thresh]
    seq_neighborhood = seq_neighborhood.replace("T", "U")
    if "U" * (u_thresh + 1) in seq_neighborhood:
        return False
    # Consecutive N
    seq_neighborhood = seq_pre[-(n_thresh):] + seq_linker + seq_post[:n_thresh]
    seq_neighborhood = seq_neighborhood.replace("T", "U")
    if any(nt * (n_thresh + 1) in seq_neighborhood for nt in set(seq_linker)):
        return False
    return True

def calc_subscores(linker_pos, *sequence_components):
    """
    Calculate base-pairing probs marginalized for each nucleotide
    """
    # Calculate bpp from ViennaRNA
    pegrna = RNA.fold_compound("".join(sequence_components))
    _ = pegrna.pf() # need to first internally calculate partition function
    basepair_probs = np.array(pegrna.bpp())[1:, 1:]
    # Fill in lower-triangle and diagonal of ViennaRNA's upper-triangular bpp matrix
    unpaired_probs = 1. - (basepair_probs.sum(axis=0) + basepair_probs.sum(axis=1))
    # copy data to make symmetric
    i_lower = np.tril_indices(len(basepair_probs), -1)
    i_diag = np.eye(len(basepair_probs), dtype=bool)
    basepair_probs[i_lower] = basepair_probs.T[i_lower]
    basepair_probs[i_diag] = unpaired_probs
    # Track indices of subsequences
    idx_cur = 0
    seq_idx = []
    for subseq in sequence_components:
        idx_prev = idx_cur
        idx_cur += len(subseq)
        seq_idx.append(slice(idx_prev, idx_cur))
    # Extract subscores for subsequences
    bpp_subseq = np.ma.masked_all(len(sequence_components))
    for i, subseq in enumerate(sequence_components):
        bpp_within_subseq = basepair_probs[seq_idx[i], seq_idx[linker_pos]]
        bpp_subseq[i] = np.mean(np.sum(bpp_within_subseq, axis=0))
    return bpp_subseq

def apply_score(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker,
                score_to_beat=None, epsilon=0.01):
    """
    Calculates subscores then outputs hashed score. Terminates calculation early if score will
    be less than score_to_beat. Prioritize PBS, spacer, template, scaffold.
    """
    # Cas9 complex at R loop subscore
    bpp_subseq1 = calc_subscores(2, seq_template, seq_pbs, seq_linker)
    subscore_pbs = 1. - bpp_subseq1[1]
    subscore_template = 1. - bpp_subseq1[0]
    # Free pegRNA subscore
    if ((score_to_beat is not None)
            and (epsilon * int(subscore_pbs / epsilon) < score_to_beat[0])):
        subscore_spacer = 0.
        subscore_scaffold = 0.
    else:
        bpp_subseq2 = calc_subscores(4, seq_spacer, seq_scaffold,
                                     seq_template, seq_pbs, seq_linker)
        subscore_spacer = 1. - bpp_subseq2[0]
        subscore_scaffold = 1. - bpp_subseq2[1]
    # Turn subscores into a single score
    return tuple(
        epsilon * int(val / epsilon)
        if val is not None else 0
        for val in (subscore_pbs, subscore_spacer, subscore_template, subscore_scaffold)
        )

def optimize(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
             linker_pattern, ac_thresh, u_thresh, n_thresh, topn, epsilon,
             num_repeats, num_steps, temp_init, temp_decay, seed):
    """
    Simulated annealing optimization of linkers
    """
    ## Pre-process inputs
    random.seed(seed)
    seq_pre = seq_spacer + seq_scaffold + seq_template + seq_pbs
    seq_post = seq_motif
    linker_pattern = linker_pattern.upper()
    ac_thresh = ac_thresh * len(linker_pattern)
    ## Simulated annealing to optimize linker sequence
    # Initialize hashmap of sequences already considered
    linker_skip = {}
    len_sequence_space = prod(len(BASE_SYMBOLS[nt]) for nt in linker_pattern)
    # Initialize min heap of topn linkers
    linker_heap = []
    for _ in range(num_repeats):
        # Initialize simulated annealing
        seq_linker_prev = "".join([random.choice(BASE_SYMBOLS[nt]) for nt in linker_pattern])
        score_prev = None
        temp = temp_init
        for _ in range(num_steps):
            # Generate new sequence by substituting characters in sequence until pass filters
            seq_linker = seq_linker_prev
            keep_going = True
            while keep_going:
                char_pos = random.randint(0, len(linker_pattern) - 1)
                seq_linker = (
                    seq_linker[:char_pos]
                    + random.choice(BASE_SYMBOLS[linker_pattern[char_pos]])
                    + seq_linker[(char_pos + 1):])
                keep_going = (
                    (seq_linker in linker_skip
                    or not apply_filters(seq_pre, seq_linker, seq_post,
                                         ac_thresh, u_thresh, n_thresh))
                    and len(linker_skip) < len_sequence_space) # already screened whole seq space
                linker_skip[seq_linker] = True
            # Calculate score for linker sequence
            score_to_beat = linker_heap[0][0] if len(linker_heap) >= topn else None
            score = apply_score(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker,
                                score_to_beat=score_to_beat, epsilon=epsilon)
            # Add to min heap i.e. maintains the top `topn` largest entries
            if score_to_beat is None: # heap is not yet full
                heapq.heappush(linker_heap, (score, seq_linker))
            elif score > score_to_beat:
                heapq.heapreplace(linker_heap, (score, seq_linker))
            # Decide if keep proposal
            if (score_prev is None                                      # initialize
                or score > score_prev                                   # exploit improvement
                or random.random() < sigmoid(                                  # explore
                    sum((s1 - s2) * (epsilon ** i)
                        for i, (s1, s2) in enumerate(zip(score, score_prev))) / temp
                    )):
                seq_linker_prev = seq_linker
                score_prev = score
            # Update simulated annealing param
            temp *= temp_decay
    linker_heap_scores, linker_heap = zip(*linker_heap)
    return linker_heap_scores, linker_heap

def apply_bottleneck(heap_scores, heap, bottleneck, seed):
    """
    Cluster sequences and output top-scoring sequence per cluster.
    """
    random.seed(seed)
    # Pick best, randomly tiebreak if needed
    def _pick_best(scores, choices):
        idx_maxed = np.where(scores == np.max(scores))[0]
        idx_chosen = random.choice(idx_maxed)
        return choices[idx_chosen]
    # Can just pick best output
    if bottleneck == 1:
        return [_pick_best(heap_scores, heap)]
    # Calculate features for each linker sequence i.e. edit distance to all other linker sequences
    features = np.zeros((len(heap), len(heap)), dtype=int)
    for i, seq_x in enumerate(heap):
        for j, seq_y in enumerate(heap):
            features[i, j] = levenshtein_distance(seq_x, seq_y)
    # Cluster linker sequences
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
    return output

def pegLIT(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
           linker_pattern="NNNNNNNN", ac_thresh=0.5, u_thresh=3, n_thresh=3, topn=100,
           epsilon=1e-2, num_repeats=10, num_steps=250, temp_init=0.15, temp_decay=0.95,
           bottleneck=1, seed=2020):
    """
    Optimizes+bottlenecks linker for an inputted pegRNA. Outputs linker recommendation(s).
    """
    # Simulated annealing to optimize linker sequence
    linker_heap_scores, linker_heap = optimize(
        seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
        linker_pattern=linker_pattern, ac_thresh=ac_thresh, u_thresh=u_thresh,
        n_thresh=n_thresh, topn=topn, epsilon=epsilon, num_repeats=num_repeats,
        num_steps=num_steps, temp_init=temp_init, temp_decay=temp_decay, seed=seed)
    # Sample diverse sequences
    linker_output = apply_bottleneck(linker_heap_scores, linker_heap,
                                     bottleneck=bottleneck, seed=seed)
    return linker_output

if __name__ == "__main__":
    # Example usage for HEK3 +1 FLAG ins
    print(pegLIT(
        seq_spacer="GGCCCAGACTGAGCACGTGA",
        seq_scaffold="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTAT"
                    "CAACTTGAAAAAGTGGCACCGAGTCGGTGC",
        seq_template="TGGAGGAAGCAGGGCTTCCTTTCCTCTGCCATCACTTATCG"
                    "TCGTCATCCTTGTAATC",
        seq_pbs="CGTGCTCAGTCTG",
        seq_motif="CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"))
