"""
Simulated annealing optimization of linkers
"""
import heapq
import random
from scipy.special import expit as sigmoid
from .score import apply_score, apply_filters
from . import constants
from .utils import sequence_space

def optimize(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
             linker_pattern=None, ac_thresh=None, u_thresh=None, n_thresh=None, topn=None,
             epsilon=None, num_repeats=None, num_steps=None, temp_init=None, temp_decay=None,
             seed=None, progress=None):
    """
    Simulated annealing optimization of linkers
    """
    ## Process arguments with defaults as needed
    random.seed(seed)
    seq_pre = seq_spacer + seq_scaffold + seq_template + seq_pbs
    seq_post = seq_motif
    if linker_pattern is None:
        linker_pattern = constants.DEFAULT_LINKER_PATTERN
    linker_pattern = linker_pattern.upper()
    if ac_thresh is None:
        ac_thresh = constants.DEFAULT_AC_THRESH
    ac_thresh = ac_thresh * len(linker_pattern)
    if u_thresh is None:
        u_thresh = constants.DEFAULT_U_THRESH
    if n_thresh is None:
        n_thresh = constants.DEFAULT_N_THRESH
    if topn is None:
        topn = constants.DEFAULT_TOPN
    if epsilon is None:
        epsilon = constants.DEFAULT_EPSILON
    if num_repeats is None:
        num_repeats = constants.DEFAULT_NUM_REPEATS
    if num_steps is None:
        num_steps = constants.DEFAULT_NUM_STEPS
    if temp_init is None:
        temp_init = constants.DEFAULT_TEMP_INIT
    if temp_decay is None:
        temp_decay = constants.DEFAULT_TEMP_DECAY
    ## Simulated annealing to optimize linker sequence
    # Initialize num linkers rejected for reasons
    SCORE_REASONS = ("PBS", "Spacer", "Template", "Scaffold")
    filter_stats = {"AC thresh": 0, "U thresh": 0, "N thresh": 0,
                    **{rsn: 0 for rsn in SCORE_REASONS}}
    # Initialize hashmap of sequences already considered
    linker_skip = {}
    len_sequence_space = sequence_space(linker_pattern)
    # Initialize min heap of topn linkers
    linker_heap = []
    for _ in range(num_repeats):
        # Initialize simulated annealing
        seq_linker_prev = "".join([random.choice(constants.BASE_SYMBOLS[nt]) for nt in linker_pattern])
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
                    + random.choice(constants.BASE_SYMBOLS[linker_pattern[char_pos]])
                    + seq_linker[(char_pos + 1):]
                    )
                if seq_linker in linker_skip:               # already screened it
                    continue
                if len(linker_skip) >= len_sequence_space:  # already screened whole seq space
                    break
                linker_skip[seq_linker] = True
                filt_pass, rsn = apply_filters(seq_pre, seq_linker, seq_post,
                                               ac_thresh, u_thresh, n_thresh, verbose=True)
                if filt_pass:
                    break
                else:
                    filter_stats[rsn] += 1
            # Calculate score for linker sequence
            score_to_beat = linker_heap[0][0] if len(linker_heap) >= topn else None
            score = apply_score(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker,
                                score_to_beat=score_to_beat, subseq_pool="mean", epsilon=epsilon)
            # Add to min heap i.e. maintains the top `topn` largest entries
            if score_to_beat is None: # heap is not yet full
                heapq.heappush(linker_heap, (score, seq_linker))
            elif score > score_to_beat:
                heapq.heapreplace(linker_heap, (score, seq_linker))
                # Track worst linker in heap is rejected
                rsn = SCORE_REASONS[[snew > sold for snew, sold in zip(score, score_to_beat)].index(True)]
                filter_stats[rsn] += 1
            else:
                # Track current linker is rejected
                rsn = SCORE_REASONS[[snew > sold for snew, sold in zip(score, score_to_beat)].index(False)]
                filter_stats[rsn] += 1
            # Decide if keep proposal
            if (score_prev is None                                          # initialize
                    or score > score_prev                                   # exploit improvement
                    or random.random() < sigmoid(                           # explore
                        sum((s1 - s2) * (epsilon ** i)
                            for i, (s1, s2) in enumerate(zip(score, score_prev))) / temp
                        )):
                seq_linker_prev = seq_linker
                score_prev = score
            # Update simulated annealing param
            temp *= temp_decay
            # Update progressbar
            if progress:
                progress.increment_step()
        if progress:
            progress.increment_repeat()
    if progress:
        progress.done()
    linker_heap_scores, linker_heap = zip(*linker_heap)
    return linker_heap_scores, linker_heap, filter_stats
