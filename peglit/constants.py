"""
Default values
"""

BASE_SYMBOLS = {
    "A": ("A",),
    "C": ("C",),
    "G": ("G",),
    "T": ("T",),
    "U": ("T",),
    "W": ("A", "T"),
    "S": ("C", "G"),
    "M": ("A", "C"),
    "K": ("G", "T"),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}

DEFAULT_LINKER_PATTERN = "NNNNNNNN"
DEFAULT_AC_THRESH = 0.5
DEFAULT_U_THRESH = 3
DEFAULT_N_THRESH = 3
DEFAULT_TOPN = 100
DEFAULT_EPSILON = 1e-2
DEFAULT_NUM_REPEATS = 10
DEFAULT_NUM_STEPS = 250
DEFAULT_TEMP_INIT = 0.15
DEFAULT_TEMP_DECAY = 0.95
DEFAULT_BOTTLENECK = 1
DEFAULT_SCAFFOLD_THRESH = 0.15
DEFAULT_MOTIF_THRESH = 0.15
DEFAULT_RNG_SEED = 2020

EXCEL_EXTENSIONS = (".xls", ".xlsx", ".xlsm", ".xlsb", ".odf", ".ods", ".odt")
