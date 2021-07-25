"""
Define scoring and functions
"""
import argparse
from .utils import make_bpp, make_bpp_subseq

def apply_filters(seq_pre, seq_linker, seq_post, ac_thresh, u_thresh, n_thresh, verbose=False):
    """
    Returns False (and filter name if verbose) if any filter is failed i.e. AC content < ac_thresh
    OR consecutive Us > u_thresh OR consecutive Ns > n_thresh. Otherwise, True if all filters are
    passed. All thresholds have units nt (i.e. ac_thresh is not a percent). Ts are treated as Us.
    """
    # AC content
    if seq_linker.count("A") + seq_linker.count("C") < ac_thresh:
        return (False, "AC thresh") if verbose else False
    # Consecutive U
    seq_neighborhood = seq_pre[-(u_thresh):] + seq_linker + seq_post[:u_thresh]
    seq_neighborhood = seq_neighborhood.replace("T", "U")
    if "U" * (u_thresh + 1) in seq_neighborhood:
        return (False, "U thresh") if verbose else False
    # Consecutive N
    seq_neighborhood = seq_pre[-(n_thresh):] + seq_linker + seq_post[:n_thresh]
    seq_neighborhood = seq_neighborhood.replace("T", "U")
    if any(nt * (n_thresh + 1) in seq_neighborhood for nt in set(seq_linker)):
        return (False, "N thresh") if verbose else False
    return (True, None) if verbose else True

def calc_subscores(linker_pos, *sequence_components, subseq_pool=None):
    """
    Calculate base-pairing probs marginalized for each nucleotide
    """
    bpp, _, seq_idx = make_bpp(*sequence_components)
    bpp_subseq = make_bpp_subseq(linker_pos, bpp, seq_idx, method=subseq_pool)
    return bpp_subseq

def apply_score(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker,
                score_to_beat=None, epsilon=0.01, subseq_pool="mean", verbose=False):
    """
    Calculates subscores and also outputs hashed score. Terminates calculation early if score will
    be less than score_to_beat. See make_bpp_subseq for info about subseq_pool. Prioritize PBS,
    spacer, template, scaffold.
    """
    # Cas9 complex at R loop
    bpp_subseq1 = calc_subscores(4, None, None, seq_template, seq_pbs, seq_linker, subseq_pool=subseq_pool)
    subscore_pbs = 1. - bpp_subseq1[3]
    subscore_template = 1. - bpp_subseq1[2]
    # Free pegRNA
    if (score_to_beat is not None) and (epsilon * int(subscore_pbs / epsilon) < score_to_beat[0]):
        subscore_spacer = None
        subscore_scaffold = None
    else:
        bpp_subseq2 = calc_subscores(4, seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker, subseq_pool=subseq_pool)
        subscore_spacer = 1. - bpp_subseq2[0]
        subscore_scaffold = 1. - bpp_subseq2[1]
    # Make scores
    score = tuple(
        epsilon * int(val / epsilon)
        if val is not None else 0
        for val in (subscore_pbs, subscore_spacer, subscore_template, subscore_scaffold)
        )
    if verbose:
        return score, (
            subscore_pbs,
            subscore_spacer,
            subscore_template,
            subscore_scaffold,
            )
    return score

def main(raw_args=None):
    parser = argparse.ArgumentParser(
        usage="%(prog)s spacer,scaffold,template,pbs linker [options]",
        description="Scores.",
        add_help=False)
    parser.add_argument("sequence", type=str)
    parser.add_argument("linker", type=str)
    parser.add_argument("--max", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(raw_args)
    print(apply_score(*args.sequence.split(","), args.linker,
                      subseq_pool=("max" if args.max else "mean"),
                      verbose=args.verbose))

if __name__ == "__main__":
    main()
