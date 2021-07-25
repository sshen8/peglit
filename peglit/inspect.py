import math
import argparse
import RNA
from .utils import make_bpp, make_bpp_subseq
from .score import apply_score
from . import constants

def make_structures(linker_pos, omit_pos, *sequence_components,
                    scaffold_pos=None, scaffold_thresh=None,
                    motif_pos=None, motif_thresh=None):
    """
    Returns predicted structures: ensemble, MFE, and MFE involving linker.
    linker_pos is index of linker sequence in sequence_components.
    """
    assert (scaffold_thresh is None) or (scaffold_pos is not None)
    assert (motif_thresh is None) or (motif_pos is not None)
    bpp, seq, seq_idx = make_bpp(*[subseq if pos not in omit_pos else None
                                   for pos, subseq in enumerate(sequence_components)])
    linker_unstructured = "." * len(seq[seq_idx[linker_pos]])
    pegrna = RNA.fold_compound(seq)
    # whole ensemble
    structure_ens, energy_ens = pegrna.pf()
    # single mfe structure
    structure_mfe, energy_mfe = pegrna.mfe()
    # suboptimal structures
    structures_subopt = pegrna.subopt_zuker()
    structure_sub, energy_sub = None, None
    for structure_subopt in structures_subopt:
        # Keep and stop if linker is structured
        if not structure_subopt.structure[seq_idx[linker_pos]] == linker_unstructured:
            structure_sub, energy_sub = structure_subopt.structure, structure_subopt.energy
            break
    # check scaffold and motif structures
    def _check_intxn(pos, thresh):
        if thresh is None or pos in omit_pos:
            return None
        bpp_subseq = make_bpp_subseq(pos, bpp, seq_idx, method="mean")
        return 1. - bpp_subseq[pos] > thresh
    scaffold_intxn = _check_intxn(scaffold_pos, scaffold_thresh)
    motif_intxn = _check_intxn(motif_pos, motif_thresh)
    # Prepare output
    def _add_filler(_seq):
        if not _seq:
            return _seq
        for pos, subseq, subseq_idx in zip(reversed(range(len(sequence_components))),
                                           reversed(sequence_components),
                                           reversed(seq_idx)):
            if pos in omit_pos:
                _seq = _seq[:subseq_idx.start] + "-" * len(subseq) + _seq[subseq_idx.stop:]
        return _seq
    return {
        "seq": _add_filler(seq),
        "scaffold_intxn": scaffold_intxn,
        "motif_intxn": motif_intxn,
        "ensemble": _add_filler(structure_ens),
        "ensemble_energy": energy_ens,
        "mfe": _add_filler(structure_mfe),
        "mfe_energy": energy_mfe,
        **({
        "subopt": _add_filler(structure_sub),
        "subopt_energy": energy_sub
        } if structure_sub else {})
    }

def calc_info(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif, seq_linker, epsilon, scaffold_thresh, motif_thresh):
    # Parse linker stats for output
    linker_info = {}
    score, subscores = apply_score(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_linker,
                                   subseq_pool="mean", epsilon=epsilon, verbose=True)
    linker_info = {
        "structures": {
            "extension": make_structures(4, (0, 1), seq_spacer, seq_scaffold, seq_template,
                                            seq_pbs, seq_linker, seq_motif, scaffold_pos=1,
                                            scaffold_thresh=scaffold_thresh,
                                            motif_pos=5, motif_thresh=motif_thresh),
            "full": make_structures(4, (), seq_spacer, seq_scaffold, seq_template, seq_pbs,
                                    seq_linker, seq_motif, scaffold_pos=1,
                                    scaffold_thresh=scaffold_thresh,
                                    motif_pos=5, motif_thresh=motif_thresh)
            },
        "score": "{:.{rnd:d}f};{:.{rnd:d}f};{:.{rnd:d}f};{:.{rnd:d}f}".format(
            *score, rnd=math.ceil(abs(math.log10(epsilon)))),
        "subscore_pbs": subscores[0],
        "subscore_spacer": subscores[1],
        "subscore_template": subscores[2],
        "subscore_scaffold": subscores[3],
        }
    return linker_info

def print_input_sequence(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif, linker_pattern):
    fmt_sequence = ("> {{spacer:=^{len_spacer}s}}"
                    "{{scaffold:-^{len_scaffold}s}}"
                    "{{template:=^{len_template}s}}"
                    "{{pbs:-^{len_pbs}s}}"
                    "{{linker:=^{len_linker}s}}"
                    "{{motif:-^{len_motif}s}}".format(len_spacer=len(seq_spacer),
                                                        len_scaffold=len(seq_scaffold),
                                                        len_template=len(seq_template),
                                                        len_pbs=len(seq_pbs),
                                                        len_linker=len(linker_pattern),
                                                        len_motif=len(seq_motif)))
    print(fmt_sequence.format(spacer="spacer", scaffold="scaffold", template="template",
                                pbs="PBS", linker="linker", motif="motif"))
    print(fmt_sequence.format(spacer=seq_spacer, scaffold=seq_scaffold, template=seq_template,
                                pbs=seq_pbs, linker=linker_pattern, motif=seq_motif))

def print_filter_stats(filter_stats, linker_pattern, ac_thresh, u_thresh, n_thresh, bottleneck, score_cutoff):
    num_remain = sum(filter_stats.values()) + bottleneck
    pad = len(str(num_remain))
    print("> [{: >{pad:}d}] All sequence space {:s}".format(num_remain, linker_pattern, pad=pad))
    FILTER_NAMES = ("Simulated annealing",
                    "AC thresh", "U thresh", "N thresh",
                    "PBS", "Spacer", "Template", "Scaffold",
                    "Bottleneck")
    FILTER_NAMES_PRINT = {
        "Simulated annealing": ("Ignored in simulated annealing", "Visited in simulated annealing"),
        "AC thresh": ("AC content < {}".format(ac_thresh), "AC content \u2265 {}".format(ac_thresh)),
        "U thresh": ("{}+ consecutive U".format(u_thresh + 1), "At most {} consecutive U".format(u_thresh)),
        "N thresh": ("{}+ consecutive N".format(n_thresh + 1), "At most {} consecutive N".format(n_thresh)),
        "PBS": ("PBS subscore \u2264 {}".format(score_cutoff[0]), "PBS subscore \u2265 {}".format(score_cutoff[0])),
        "Spacer": ("Spacer subscore \u2264 {}".format(score_cutoff[1]), "Spacer subscore \u2265 {}".format(score_cutoff[1])),
        "Template": ("Template subscore \u2264 {}".format(score_cutoff[2]), "Template subscore \u2265 {}".format(score_cutoff[2])),
        "Scaffold": ("Scaffold subscore \u2264 {}".format(score_cutoff[3]), "Scaffold subscore \u2265 {}".format(score_cutoff[3])),
        "Bottleneck": ("Bottleneck", "Final output"),
    }
    for i, filt_name in enumerate(FILTER_NAMES):
        num_remain -= filter_stats[filt_name]
        print("> " + " "*((4*i)+2) + "\u251C\u2500[{num: >{pad:}d}] {title:s}".format(title=FILTER_NAMES_PRINT[filt_name][0], num=filter_stats[filt_name], pad=pad))
        print("> " + " "*((4*i)+2) + "\u2514\u2500[{num: >{pad:}d}] {title:s}".format(title=FILTER_NAMES_PRINT[filt_name][1], num=num_remain, pad=pad))

def print_structures(linker_stat, verbose):
    # Print warnings about scaffold or motif interactions
    if any(struct["scaffold_intxn"] for struct in linker_stat["structures"].values()):
        print("* Significant scaffold:pegRNA interaction predicted")
    if any(struct["motif_intxn"] for struct in linker_stat["structures"].values()):
        print("+ Significant motif:pegRNA interaction predicted")
    # Print structures and scores
    if verbose:
        fmt_scores = ("> {full_sequence:s} (score = {score:s})")
        fmt_structure = ("> {{{fold_name}:s}} "
                            "(\N{GREEK CAPITAL LETTER DELTA}G = {{{fold_name}_energy:+4.1f}} kcal/mol)"
                            )
        FOLD_NAMES = ("ensemble", "mfe", "subopt")
        fmt_structures = "\n".join([fmt_structure.format(fold_name=fold_name)
                                    for fold_name in FOLD_NAMES])
        for _, structure in linker_stat["structures"].items():
            print(fmt_scores.format(full_sequence=structure["seq"],
                                    score=linker_stat["score"]))
            print(fmt_structures.format(**structure))

def main(raw_args=None):
    parser = argparse.ArgumentParser(
        usage="%(prog)s spacer,scaffold,template,pbs,motif linker [options]",
        description="Scores.",
        add_help=False)
    parser.add_argument("sequence", type=str)
    parser.add_argument("linker", type=str)
    parser.add_argument("--epsilon", type=float, default=constants.DEFAULT_EPSILON)
    parser.add_argument("--scaffold-thresh", type=float, default=constants.DEFAULT_SCAFFOLD_THRESH)
    parser.add_argument("--motif-thresh", type=float, default=constants.DEFAULT_MOTIF_THRESH)
    args = parser.parse_args(raw_args)
    seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif = args.sequence.split(",")
    linker_stat = calc_info(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif, args.linker,
                            epsilon=args.epsilon, scaffold_thresh=args.scaffold_thresh, motif_thresh=args.motif_thresh)
    # Print recommended linker
    print(args.linker)
    # Print warnings and scores and structures if requested
    print_structures(linker_stat, verbose=True)

if __name__ == "__main__":
    main()
