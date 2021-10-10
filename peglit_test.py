import random
from peglit.score import apply_filters as filt_main
from peglit_min import apply_filters as filt_min
from peglit.score import apply_score as score_main
from peglit_min import apply_score as score_min
from peglit.bottleneck import apply_bottleneck as bottleneck_main
from peglit_min import apply_bottleneck as bottleneck_min
from peglit.optimize import optimize as optimize_main
from peglit_min import optimize as optimize_min
from peglit.peglit import pegLIT as pegLIT_main
from peglit_min import pegLIT as pegLIT_min
from peglit.utils import TqdmProgressObserver
import pytest

@pytest.mark.parametrize("case_in,case_out", [
    # AC thresh
    ({"seq_pre": "", "seq_linker": "AGGAAGGA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "", "seq_linker": "AGGCCGGA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "", "seq_linker": "AGGAGGGA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "AC thresh")),
    ({"seq_pre": "", "seq_linker": "AGGAGGGA", "seq_post": "", "ac_thresh": 3, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    # U thresh
    ({"seq_pre": "", "seq_linker": "AATATTTA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "", "seq_linker": "AATTTTAA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "U thresh")),
    ({"seq_pre": "", "seq_linker": "AATTUTAA", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "U thresh")),
    ({"seq_pre": "", "seq_linker": "AATTUTAA", "seq_post": "", "ac_thresh": 4, "u_thresh": 4, "n_thresh": 3}, (False, "N thresh")),
    ({"seq_pre": "", "seq_linker": "AATTUTAA", "seq_post": "", "ac_thresh": 4, "u_thresh": 4, "n_thresh": 4}, (True, None)),
    # N thresh
    ({"seq_pre": "", "seq_linker": "AAAGAATG", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "", "seq_linker": "AAAAGATG", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "N thresh")),
    ({"seq_pre": "", "seq_linker": "AAAAGATG", "seq_post": "", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 4}, (True, None)),
    # Context for consecutive nt filters
    ({"seq_pre": "A", "seq_linker": "TAAATAAA", "seq_post": "G", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "G", "seq_linker": "TAAATAAA", "seq_post": "G", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (True, None)),
    ({"seq_pre": "A", "seq_linker": "AAATTAAA", "seq_post": "G", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "N thresh")),
    ({"seq_pre": "G", "seq_linker": "AAATTAAA", "seq_post": "A", "ac_thresh": 4, "u_thresh": 3, "n_thresh": 3}, (False, "N thresh")),
])
def test_filters(case_in, case_out):
    assert filt_main(**case_in, verbose=True) == case_out
    assert filt_main(**case_in, verbose=False) == case_out[0]
    assert filt_min(**case_in) == case_out[0]

@pytest.mark.parametrize("case_in,case_out", [
    ({
        "seq_spacer": "GGCCCAGACTGAGCACGTGA",
        "seq_scaffold": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
        "seq_template": "TGGAGGAAGCAGGGCTTCCTTTCCTCTGCCATCACTTATCGTCGTCATCCTTGTAATC",
        "seq_pbs": "CGTGCTCAGTCTG",
        "seq_linker": "AA",
    }, (0.97, 0.99, 0.75, 0.99)),
    ({
        "seq_spacer": "GGCCCAGACTGAGCACGTGA",
        "seq_scaffold": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
        "seq_template": "TGGAGGAAGCAGGGCTTCCTTTCCTCTGCCATCACTTATCGTCGTCATCCTTGTAATC",
        "seq_pbs": "CGTGCTCAGTCTG",
        "seq_linker": "ATCGATCG",
    }, (0.99, 0.93, 0.34, 0.99)),
])
def test_scores(case_in, case_out):
    assert score_main(**case_in, verbose=True)[0] == case_out
    assert score_main(**case_in, verbose=False) == case_out
    assert score_min(**case_in) == case_out
    # motif doesn't affect score
    with pytest.raises(Exception):
        score_main(**case_in, seq_motif="ATGC")
    with pytest.raises(Exception):
        score_min(**case_in, seq_motif="ATGC")

random.seed(2020)
@pytest.mark.parametrize("case_in,case_out", (
    ({
        "heap_scores": [(1,), (0,), (1,), (0,), (1,)],
        "heap": ["AAAA", "AAAT", "TTGC", "TTGG", "TTCG"],
        "bottleneck": 2,
        "seed": 2020,
    }, {"AAAA", "TTGC"}),
    ({ # different seed to test random tie-breaking
        "heap_scores": [(1,), (0,), (1,), (0,), (1,)],
        "heap": ["AAAA", "AAAT", "TTGC", "TTGG", "TTCG"],
        "bottleneck": 2,
        "seed": 2021,
    }, {"AAAA", "TTCG"}),
    ({ # special handling for bottleneck=1
        "heap_scores": [(1,),] * 100,
        "heap": ["".join([random.choice(("A", "T", "G", "C")) for _ in range(4)]) for _ in range(100)],
        "bottleneck": 1,
        "seed": 2020,
    }, {"CTTG",}),
))
def test_bottleneck(case_in, case_out):
    assert len(set(bottleneck_main(**case_in, verbose=True)[0]) & case_out) == len(case_out)
    assert len(set(bottleneck_main(**case_in, verbose=False)) & case_out) == len(case_out)
    assert len(set(bottleneck_min(**case_in)) & case_out) == len(case_out)

@pytest.mark.parametrize("case_in,case_out", (
    ({
        "seq_spacer": "GGCCCAGACTGAGCACGTGA",
        "seq_scaffold": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
        "seq_template": "TGGAGGAAGCAGGGCTTCCTTTCCTCTGCCATCACTTATCGTCGTCATCCTTGTAATC",
        "seq_pbs": "CGTGCTCAGTCTG",
        "seq_motif": "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA",
        "linker_pattern": "NNNNnNNN", # test lowercase in pattern
        "ac_thresh": 0.5,
        "u_thresh": 3,
        "n_thresh": 3,
        "topn": 5,
        "epsilon": 0.01,
        "num_repeats": 3,
        "num_steps": 10,
        "temp_init": 0.15,
        "temp_decay": 0.95,
        "seed": 2020,
    }, {"ACTCTATG", "ACTATGAT", "ACTACGAG", "ACTATGAG", "ACTAAGAT"}),
))
def test_optimize(case_in, case_out):
    assert len(set(optimize_main(**case_in)[1]) & case_out) == len(case_out)
    assert len(set(optimize_min(**case_in)[1]) & case_out) == len(case_out)

@pytest.mark.parametrize("case_in,case_out", (
    ({ # DNMT1 +1 flag insertion
        "seq_spacer": "GATTCCTGGTGCCAGAAACA",
        "seq_scaffold": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
        "seq_template": "TCTGCCCTCCCGTCACCCCTGT",
        "seq_pbs": "TTCTGGCACCAGGA",
        "seq_motif": "GGGTCAGGAGCCCCCCCCCTGAACCCAGGATAACCCTCAAAGTCGGGGGGCAACCC",
        "linker_pattern": "NNNNNNNNNNNNNNNNNN",
        "topn": 10,
        "num_repeats": 3,
        "num_steps": 5,
        "seed": 2020,
    }, {"ACCGTAAACATTAAGGTT",}),
))
def test_pegLIT(case_in, case_out):
    assert len(set(pegLIT_main(**case_in, verbose=False)) & case_out) == len(case_out)
    assert len(pegLIT_main(**case_in, verbose=True)) == 5
    assert len(set(pegLIT_min(**case_in)) & case_out) == len(case_out)
    assert len(set(pegLIT_main(**case_in, progress=TqdmProgressObserver(case_in["num_repeats"], case_in["num_steps"]))) & case_out) == len(case_out)
