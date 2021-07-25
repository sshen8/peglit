"""
Main entry point for running optimization+bottlenecking. Accepts either an individual or a CSV
batch of pegRNA sequences; output is printed or written to a new CSV, respectively.
"""
import os
import argparse
from time import sleep
import pandas as pd
from tqdm import tqdm
from matplotlib import pyplot as plt
from .utils import sequence_space, TqdmProgressObserver
from .plots import plot_clusters
from .optimize import optimize
from .bottleneck import apply_bottleneck
from .inspect import calc_info, print_input_sequence, print_filter_stats, print_structures
from . import constants

def pegLIT(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
           linker_pattern="NNNNNNNN", ac_thresh=0.5, u_thresh=3, n_thresh=3, topn=100,
           epsilon=1e-2, num_repeats=10, num_steps=250, temp_init=0.15, temp_decay=0.95,
           bottleneck=1, seed=2020, verbose=False, progress=None):
    """
    Optimizes+bottlenecks linker for an inputted pegRNA. Outputs linker recommendation(s).
    """
    # Simulated annealing to optimize linker sequence
    linker_heap_scores, linker_heap, filter_stats = optimize(
        seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
        linker_pattern=linker_pattern, ac_thresh=ac_thresh, u_thresh=u_thresh,
        n_thresh=n_thresh, topn=topn, epsilon=epsilon, num_repeats=num_repeats,
        num_steps=num_steps, temp_init=temp_init, temp_decay=temp_decay, seed=seed,
        progress=progress)
    # Sample diverse sequences
    linker_output, linker_feats = apply_bottleneck(linker_heap_scores, linker_heap,
                                                   bottleneck=bottleneck, seed=seed, verbose=True)
    if verbose:
        return linker_output, linker_feats, linker_heap_scores, linker_heap, filter_stats
    return linker_output

def make_output(args, seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif):
    """
    Takes parsed inputs and does optimization+bottlenecking for an individual pegRNA. Outputs a
    DataFrame, each row is a recommended linker.
    """
    ## Do pegLIT
    progress = TqdmProgressObserver(args.num_repeats, args.num_steps)
    linker_output, linker_feats, linker_heap_scores, linker_heap, filter_stats = pegLIT(
        seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
        linker_pattern=args.linker_pattern, ac_thresh=args.ac_thresh, u_thresh=args.u_thresh,
        n_thresh=args.n_thresh, topn=args.topn, epsilon=args.epsilon, num_repeats=args.num_repeats,
        num_steps=args.num_steps, temp_init=args.temp_init, temp_decay=args.temp_decay,
        bottleneck=args.bottleneck, seed=args.seed, verbose=True, progress=progress)
    filter_stats["Bottleneck"] = args.topn - args.bottleneck
    filter_stats["Simulated annealing"] = (sequence_space(args.linker_pattern)
                                           - sum(filter_stats.values())
                                           - args.bottleneck)
    ## Make outputs
    # Make plots of clusters if requested
    if args.clusters is not None:
        plot_clusters(linker_heap, linker_feats, args.bottleneck)
        plt.savefig(args.clusters.name, bbox_inches="tight")
    return linker_output, filter_stats, linker_heap_scores[0]

def main(raw_args=None):
    """
    Main entry point for running optimization+bottlenecking.
    """
    # Process inputs
    parser = argparse.ArgumentParser(
        usage="%(prog)s (sequence | --batch FILE) [options]",
        description="Searches for linker sequences, avoiding complementarity to pegRNA sequence.",
        add_help=False)
    group_seq = parser.add_argument_group(title="pegRNA sequence(s)").add_mutually_exclusive_group(required=True)
    group_seq.add_argument("sequence", nargs='?',
                           help="A single pegRNA sequence, to which linker should avoid "
                                "complementarity, in 5' -> 3' orientation. Comma-separated "
                                "subsequences denote: spacer, scaffold, template, PBS, motif.")
    group_seq.add_argument("--batch", metavar="FILE", type=argparse.FileType("r+"),
                           help="A headered CSV file containing pegRNA sequences. Must contain "
                                "columns for spacer, scaffold, template, PBS, motif. Output will be "
                                "written in a CSV with the same name but appended with "
                                "_linker_designs.")
    group_par = parser.add_argument_group(title="algorithm parameters")
    group_par.add_argument("--linker-pattern", type=str, default=constants.DEFAULT_LINKER_PATTERN,
                           help="Nucleotides allowed in linker. Default: {}."
                           .format(constants.DEFAULT_LINKER_PATTERN))
    group_par.add_argument("--ac-thresh", type=float, default=constants.DEFAULT_AC_THRESH,
                           help="Minimum allowed AC content in linker sequence. Default: {}."
                           .format(constants.DEFAULT_AC_THRESH))
    group_par.add_argument("--u-thresh", type=int, default=constants.DEFAULT_U_THRESH,
                           help="Maximum number of consecutive Us allowed in linker sequence."
                                "Default: {}."
                           .format(constants.DEFAULT_U_THRESH))
    group_par.add_argument("--n-thresh", type=int, default=constants.DEFAULT_N_THRESH,
                           help="Maximum number of consecutive any nucleotide allowed."
                                "Default: {}."
                           .format(constants.DEFAULT_N_THRESH))
    group_par.add_argument("--topn", type=int, default=constants.DEFAULT_TOPN,
                           help="Keep this many of the best linkers. Small value -> better "
                                "linker sequences. Large value -> potentially more diverse. "
                                "Default: {}."
                           .format(constants.DEFAULT_TOPN))
    group_par.add_argument("--epsilon", type=float, default=constants.DEFAULT_EPSILON,
                           help="Basepairing probabilities are considered equal if their "
                                "difference is less than this tolerance value. Default: {}."
                           .format(constants.DEFAULT_EPSILON))
    group_par.add_argument("--seed", type=int, default=constants.DEFAULT_RNG_SEED,
                           help="Reproducibly initializes pseudorandom number generator. "
                                "Default: {}."
                           .format(constants.DEFAULT_RNG_SEED))
    group_anl = parser.add_argument_group(title="simulated annealing parameters")
    group_anl.add_argument("--num-repeats", type=int, default=constants.DEFAULT_NUM_REPEATS,
                           help="Number of repeats. Default: {}."
                           .format(constants.DEFAULT_NUM_REPEATS))
    group_anl.add_argument("--num-steps", type=int, default=constants.DEFAULT_NUM_STEPS,
                           help="Number of steps. Default: {}."
                           .format(constants.DEFAULT_NUM_STEPS))
    group_anl.add_argument("--temp-init", type=float, default=constants.DEFAULT_TEMP_INIT,
                           help="Initial temperature. Default: {}."
                           .format(constants.DEFAULT_TEMP_INIT))
    group_anl.add_argument("--temp-decay", type=float, default=constants.DEFAULT_TEMP_DECAY,
                           help="Multiplicative temperature decay per step. Default: {}."
                           .format(constants.DEFAULT_TEMP_DECAY))
    group_out = parser.add_argument_group(title="output options")
    group_out.add_argument("--bottleneck", type=int, default=constants.DEFAULT_BOTTLENECK,
                           help="Number of sequences to output. Default: {}."
                           .format(constants.DEFAULT_BOTTLENECK))
    group_out.add_argument("--verbose", action="store_true",
                           help="In addition to linker sequences, print MFE structure and top "
                                "Zucker suboptimal structure with linker complementarity.")
    group_out.add_argument("--scaffold-thresh", type=float, default=constants.DEFAULT_SCAFFOLD_THRESH,
                           help="In verbose output, * indicates predicted structure contains "
                                "scaffold that interacts with the rest of the pegRNA with "
                                "probability greater than threshold. Default: {}."
                           .format(constants.DEFAULT_SCAFFOLD_THRESH))
    group_out.add_argument("--motif-thresh", type=float, default=constants.DEFAULT_MOTIF_THRESH,
                           help="In verbose output, + indicates predicted structure contains "
                                "motif that interacts with the rest of the pegRNA with "
                                "probability greater than threshold. Default: {}."
                           .format(constants.DEFAULT_SCAFFOLD_THRESH))
    group_out.add_argument("--clusters", type=argparse.FileType("w"),
                           help="Filename to save clusters plot.")
    group_out.add_argument("-h", "--help", action="help",
                           help="Show this help message and exit")
    group_col = parser.add_argument_group(title="(advanced) set algorithm parameters per pegRNA in CSV")
    group_col.add_argument("--linker-pattern-col", type=str)
    group_col.add_argument("--ac-thresh-col", type=str)
    group_col.add_argument("--u-thresh-col", type=str)
    group_col.add_argument("--n-thresh-col", type=str)
    group_col.add_argument("--topn-col", type=str)
    group_col.add_argument("--epsilon-col", type=str)
    group_col.add_argument("--num-repeats-col", type=str)
    group_col.add_argument("--num-steps-col", type=str)
    group_col.add_argument("--temp-init-col", type=str)
    group_col.add_argument("--temp-decay-col", type=str)
    group_col.add_argument("--bottleneck-col", type=str)
    group_col.add_argument("--scaffold-thresh-col", type=str)
    group_col.add_argument("--motif-thresh-col", type=str)
    args = parser.parse_args(raw_args)
    assert 0 < args.bottleneck < args.topn, "Bottleneck needs to be smaller than topn."
    # Individual input. Output is printed
    if args.sequence is not None:
        seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif = args.sequence.split(",")
        pbar = tqdm(total=1, position=0, desc="pegRNA")
        linker_output, filter_stats, score_cutoff = make_output(args, seq_spacer, seq_scaffold,
                                                                seq_template, seq_pbs, seq_motif)
        linker_info = {seq_linker: calc_info(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif, seq_linker,
                                             epsilon=args.epsilon, scaffold_thresh=args.scaffold_thresh, motif_thresh=args.motif_thresh) for seq_linker in linker_output}
        pbar.update()
        pbar.close()
        # Print inputted info
        print_input_sequence(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif, args.linker_pattern)
        # Print outputted info
        if args.verbose:
            print_filter_stats(filter_stats, args.linker_pattern, args.ac_thresh, args.u_thresh, args.n_thresh, args.bottleneck, score_cutoff)
        for seq_linker, linker_stat in linker_info.items():
            # Print recommended linker
            print(seq_linker)
            # Print warnings and scores and structures if requested
            print_structures(linker_stat, verbose=args.verbose)
    # Batch input read from CSV file. Output into new columns in duplicated CSV file
    elif args.batch is not None:
        # Make and validate output filename
        fname_name, fname_xt = os.path.splitext(args.batch.name)
        fname_xt = fname_xt.lower()
        supported_xt = constants.EXCEL_EXTENSIONS + (".csv",)
        if fname_xt not in supported_xt:
            raise AssertionError("Please make sure your batch input is an Excel or CSV file. "
                                 "Your extension: {:s}. Supported extensions: {:s}.".format(
                                     fname_xt, ", ".join(supported_xt)
                                 ))
        fname_output = "{:s}_linker_designs{:s}".format(fname_name, fname_xt)
        index = 2
        while os.path.exists(fname_output):
            fname_output = "{:s}_linker_designs_{:d}{:s}".format(fname_name, index, fname_xt)
            index += 1
        # Read and process inputs
        err = []
        for reader in (pd.read_excel, pd.read_csv):
            try:
                df_input = reader(args.batch.name).dropna(how="all").fillna("")
                break
            except Exception as _err:
                df_input = None
                err.append(_err)
                continue
        if df_input is None:
            raise AssertionError("Could not read file as Excel or CSV: {:s}.\n{:s}".format(args.batch.name, "\n".join(_err.__str__() for _err in err)))
        if len(df_input) > 40:
            print("Your batch of {:d} pegRNAs might take a while. Make sure the session will stay "
                  "alive long enough. Linker generation starts in 2 sec.".format(len(df_input)))
            sleep(2)
        df_output = []
        for _, sequence in tqdm(df_input.iterrows(), position=0, total=len(df_input), desc="pegRNA"):
            # Get user parameters for individual pegRNA (overrides general parameters)
            if args.linker_pattern_col is not None:
                args.linker_pattern = sequence[args.linker_pattern_col]
            if args.ac_thresh_col is not None:
                args.ac_thresh = sequence[args.ac_thresh_col]
            if args.u_thresh_col is not None:
                args.u_thresh = sequence[args.u_thresh_col]
            if args.n_thresh_col is not None:
                args.n_thresh = sequence[args.n_thresh_col]
            if args.topn_col is not None:
                args.topn = sequence[args.topn_col]
            if args.epsilon_col is not None:
                args.epsilon = sequence[args.epsilon_col]
            if args.num_repeats_col is not None:
                args.num_repeats = sequence[args.num_repeats_col]
            if args.num_steps_col is not None:
                args.num_steps = sequence[args.num_steps_col]
            if args.temp_init_col is not None:
                args.temp_init = sequence[args.temp_init_col]
            if args.temp_decay_col is not None:
                args.temp_decay = sequence[args.temp_decay_col]
            if args.bottleneck_col is not None:
                args.bottleneck = sequence[args.bottleneck_col]
            if args.scaffold_thresh_col is not None:
                args.scaffold_thresh = sequence[args.scaffold_thresh_col]
            if args.motif_thresh_col is not None:
                args.motif_thresh = sequence[args.motif_thresh_col]
            # Run everything for individual pegRNA
            linker_output, filter_stats, score_cutoff = make_output(args, sequence["spacer"],
                                                                    sequence["scaffold"],
                                                                    sequence["template"],
                                                                    sequence["PBS"],
                                                                    sequence["motif"])
            linker_info = {
                seq_linker: calc_info(sequence["spacer"], sequence["scaffold"],
                                      sequence["template"], sequence["PBS"], sequence["motif"],
                                      seq_linker, epsilon=args.epsilon,
                                      scaffold_thresh=args.scaffold_thresh,
                                      motif_thresh=args.motif_thresh)
                for seq_linker in linker_output}
            for seq_linker, linker_stat in linker_info.items():
                sequence_new = sequence.copy(deep=True)
                sequence_new["linker"] = seq_linker
                if args.verbose:
                    sequence_new = sequence_new.append(
                        pd.json_normalize(linker_stat, sep="_").squeeze())
                    sequence_new = sequence_new.append(
                        pd.json_normalize({"filter_{}".format(name): num for name, num in filter_stats.items()}, sep="_").squeeze())
                    sequence_new["filter_scorecutoff"] = ";".join(str(subscore) for subscore in score_cutoff)
                df_output.append(sequence_new)
        # Write output to new file
        df_output = pd.DataFrame(df_output)
        if fname_xt in constants.EXCEL_EXTENSIONS:
            df_output.to_excel(fname_output, index=False)
        elif fname_xt == ".csv":
            df_output.to_csv(fname_output, index=False)
        else:
            raise AssertionError()
        print(fname_output)
    else:
        raise ValueError()

if __name__ == "__main__":
    main()
