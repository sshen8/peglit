# pegRNA Linker Identification Tool (pegLIT)

pegRNA Linker Identification Tool (pegLIT) automatically identifies non-interfering nucleotide linkers between a pegRNA and 3' motif. For more about the underlying science, please check out [our publication](https://doi.org/10.1038/s41587-021-01039-7):

    J. Nelson, P. Randolph, S. Shen, K. Everette, P. Chen, A. Anzalone, M. An, G. Newby, J. Chen, A. Hsu, and D. Liu. Engineered pegRNAs that improve prime editing efficiency. Nature Biotechnology (2021).

## Installation

There are two main ways to access and use pegLIT.

* The **web app** is probably the easiest and is sufficiently powerful for most use cases. You can access it at [peglit.liugroup.us](https://peglit.liugroup.us).
* If your pegRNAs are unusually long (say >150 nt), pegLIT might need more time than the 10-minute limit we've imposed for our server. In that case, you can install pegLIT onto your computer via our **Python package**, which will let you programatically run pegLIT as part of your own custom pipeline.

### Python package

To install pegLIT from [Bioconda](https://bioconda.github.io) into your current conda environment:
```
conda install peglit
```

#### Additional steps for Apple silicon
If you're running a [Mac with Apple silicon](https://support.apple.com/en-us/HT211814), you need to do a few more things because some of pegLIT's dependencies aren't compatible with Apple silicon yet.
1. Install [Rosetta](https://support.apple.com/en-us/HT211861) if you haven't already:
    ```
    softwareupdate --install-rosetta
    ```
2. Create and activate a conda environment that uses the Intel versions of packages:
    ```
    CONDA_SUBDIR=osx-64 conda create -n peglit_env
    conda activate peglit_env
    ```
3. Now you can install pegLIT.

## Usage

### Command Line Interface
```
peglit ATGC,ATGC,ATGC,ATGC,ATGC
```
The input is a pegRNA sequence with commas separating the spacer, scaffold, template, PBS, and motif subsequences, respectively. The output will be printed.

You can also run a _batch_ of pegRNA sequences in a CSV file:
```
peglit --batch batch_pegRNA_sequences.csv
```
The first row of the CSV file should include the headers: `spacer`, `scaffold`, `template`, `PBS`, `motif`; and each subsequent row should contain the corresponding pegRNA subsequences. The output will be saved to a new CSV file: `<input_filename>_linker_designs.csv`.

### Python

```
import peglit

linkers = peglit.pegLIT(seq_spacer="ATGC", seq_scaffold="ATGC", seq_template="ATGC",
                        seq_pbs="ATGC", seq_motif="ATGC")
```
The output `linkers` is a list of recommended linker sequence(s). In this case, `['TCGACCCT']`.

To calculate linker basepairing scores for an epegRNA:
```
import peglit

subscores = peglit.score(seq_spacer="ATGC", seq_scaffold="ATGC", seq_template="ATGC",
                         seq_pbs="ATGC", seq_linker="ATGC")
```
The output `subscores` is a tuple containing the PBS, spacer, template, and scaffold subscores. In this case, `(0.90, 0.51, 0.69, 0.55)`. Each subscore is between 0 and 1; higher is better (less basepairing).

## Documentation

```
def pegLIT(seq_spacer, seq_scaffold, seq_template, seq_pbs, seq_motif,
           linker_pattern="NNNNNNNN", ac_thresh=0.5, u_thresh=3, n_thresh=3, topn=100,
           epsilon=1e-2, num_repeats=10, num_steps=250, temp_init=0.15, temp_decay=0.95,
           bottleneck=1, seed=2020, verbose=False, progress=None):
    """
    Optimizes+bottlenecks linker for an inputted pegRNA. Outputs linker recommendation(s).
    """
```

### Input

* `seq_spacer`, `seq_scaffold`, `seq_pbs`, `seq_motif` : string
    * Sequences of the epegRNA's spacer, scaffold, PBS, and 3' motif, respectively, in 5'-to-3' orientation. T and U are interchangeable. Case insensitive.
* `linker_pattern` : string, optional
    * Pattern that must be matched by the outputted linker sequence. Use IUPAC degenerate base symbols. Case insensitive.
    * Default: NNNNNNNN
* `ac_thresh` : float, optional
    * Minimum number of A or C bases allowed in the outputted linker, expressed as a fraction of the length of the linker.
    * Default: 0.5
* `u_thresh` : int, optional
    * Maximum number of consecutive U bases allowed in the epegRNA containing the outputted linker.
    * Default: 3
* `n_thresh` : int, optional
    * Maximum number of consecutive bases of any base allowed in the epegRNA containing the outputted linker.
    * Default: 3
* `topn` : int, optional
    * Number of "hits" to keep after simulated annealing and before bottlenecking. A larger value may increase diversity of the outputted linkers but may be more computationally expensive (because the heap of "hits" would be larger during simulated annealing, and there'd be more linkers to cluster during bottlenecking).
    * Default: 100
* `epsilon` : float, optional
    * The algorithm rounds linker basepairing subscores to the nearest `epsilon`. A larger value would cause the algorithm to pay more attention to epegRNA components that are earlier in the priority order (PBS > spacer > template > scaffold). Another perspective is that `epsilon` is the smallest meaningful difference between subscores — subtler differences are treated as prediction error.
    * Default: 0.01
* `num_repeats` : int, optional
    * Number of random re-initializations of simulated annealing. A larger value might give the algorithm more chances to find and output higher-scoring linkers at the cost of more computation time.
    * Default: 10
* `num_steps` : int, optional
    * Number of steps taken during each simulated annealing run. The tradeoff for this parameter is the same as for `num_repeats`.
    * Default: 250
* `temp_init` : float, optional
    * Initial temperature for each simulated annealing run. A larger value encourages global exploration over local optimization.
    * Default: 0.15
* `temp_decay` : float, optional
    * Temperature decay factor applied after each simulated annealing step. $Temperature(t) = (temp\_init) * (temp\_decay)^t$.
    * Default: 0.95
* `bottleneck` : int, optional
    * Number of linkers to sample and output among the "hits" from simulated annealing.
    * Default: 1
* `seed` : int, optional
    * Seed to random number generator. Good to pick a known, fixed value for reproducibility.
    * Default: 2020
* `verbose` : bool, optional
    * Additional outputs. See below.
    * Default: False
* `progress` : `peglit.utils.ProgressObserver`, optional
    * Progress tracking and reporting.
    * Default: None

### Output

* `linkers` : list of strings
    * List containing `bottleneck` recommended linker sequences.
* (if `verbose`) `linker_feats` : Numpy array with shape (`topn`, `topn`)
    * Pairwise Levenshtein distances between linker sequence "hits" from simulated annealing.
* (if `verbose`) `linker_heap_scores` : list of (float, float, float, float)
    * Score tuples for all linker sequence "hits" from simulated annealing.
* (if `verbose`) `linker_heap` : list of strings
    * Linker sequence "hits" from simulated annealing.
* (if `verbose`) `filter_stats` : dict
    * Keys indicate the stage at which a candidate linker was removed from consideration. Values indicate the number of such linkers.

##  Contact

We appreciate all questions, bug reports, feature requests, and other feedback about pegLIT. Please [email us](https://peglit.liugroup.us/about).