"""
Helper functions that interface with ViennaRNA.
"""
import numpy as np
from tqdm import tqdm
import RNA
from . import constants

def make_sequence(*sequence_components):
    """
    Concatenates optional subsequences. Also returns indices of subsequences.
    """
    idx_cur = 0
    seq = ""
    seq_idx = []
    for subseq in sequence_components:
        if subseq is not None:
            idx_prev = idx_cur
            idx_cur += len(subseq)
            seq += subseq
            seq_idx.append(slice(idx_prev, idx_cur))
        else:
            seq_idx.append(slice(idx_cur, idx_cur))
    return seq, seq_idx

def make_bpp(*sequence_components):
    """
    Returns matrix of base pairing probabilities from ViennaRNA calculation.
    Matrix is upper-triangular.
    """
    seq, seq_idx = make_sequence(*sequence_components)
    pegrna = RNA.fold_compound(seq)
    _ = pegrna.pf() # need to first internally calculate partition function
    basepair_probs = np.array(pegrna.bpp())[1:, 1:]
    return basepair_probs, seq, seq_idx

def make_unpaired_probs(basepair_probs):
    """
    Returns probability that each nucleotide is unpaired, using base pairing probabilities from
    ViennaRNA calculation.
    """
    unpaired_probs = 1. - (basepair_probs.sum(axis=0) + basepair_probs.sum(axis=1))
    return unpaired_probs

def fill_bpp(basepair_probs):
    """
    Fills in lower-triangle and diagonal of ViennaRNA's upper-triangular bpp matrix.
    """
    unpaired_probs = make_unpaired_probs(basepair_probs)
    # copy data to make symmetric
    i_lower = np.tril_indices(len(basepair_probs), -1)
    i_diag = np.eye(len(basepair_probs), dtype=bool)
    bpp = basepair_probs
    bpp[i_lower] = basepair_probs.T[i_lower]
    bpp[i_diag] = unpaired_probs
    return bpp

def make_bpp_subseq(linker_pos, basepair_probs, seq_idx, method="mean"):
    """
    Pools bpp info from the level of individual basepairs to the level of subsequences.
    """
    bpp_full = fill_bpp(basepair_probs.copy())
    bpp_subseq = np.ma.masked_all(len(seq_idx))
    for i, _ in enumerate(seq_idx):
        if seq_idx[i].stop - seq_idx[i].start == 0:
            continue
        bpp_within_subseq = bpp_full[seq_idx[i], seq_idx[linker_pos]]
        if method == "mean":
            bpp_subseq[i] = np.mean(np.sum(bpp_within_subseq, axis=0))
        elif method == "max":
            bpp_subseq[i] = np.max(np.sum(bpp_within_subseq, axis=0))
        else:
            raise NotImplementedError()
    return bpp_subseq

def sequence_space(pattern):
    """
    Number of sequences consistent with the pattern
    """
    return np.prod([len(constants.BASE_SYMBOLS[nt]) for nt in pattern])

class ProgressObserver:
    def __init__(self, num_repeats, num_steps):
        self.repeat_i = 0
        self.repeat_total = num_repeats
        self.step_i = 0
        self.step_total = num_steps

    def increment_step(self):
        self.step_i += 1
    
    def increment_repeat(self):
        self.repeat_i += 1
        self.step_i = 0
    
    def done(self):
        self.step_i = self.step_total
        self.repeat_i = self.repeat_total

class TqdmProgressObserver(ProgressObserver):
    def __init__(self, num_repeats, num_steps):
        super().__init__(num_repeats, num_steps)
        self.repeat_pbar = tqdm(total=self.repeat_total, position=1, desc="Repeats", leave=False)
        self.step_pbar = tqdm(total=self.step_total, position=2, desc="Steps", leave=False)
    
    def increment_step(self):
        out = super().increment_step()
        self.step_pbar.update()
        return out
    
    def increment_repeat(self):
        out = super().increment_repeat()
        self.repeat_pbar.update()
        self.step_pbar.reset()
        return out

    def done(self):
        out = super().done()
        self.repeat_pbar.close()
        self.step_pbar.close()
        return out
