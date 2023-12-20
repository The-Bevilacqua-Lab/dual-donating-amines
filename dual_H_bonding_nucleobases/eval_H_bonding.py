"""
This module will eventaully do something.
"""

from pymol import cmd
from pymol import stored


# adpated from pdb_miner.py from rare-nucleobase-miner repository
residues = {
    "ARG": {
        "donors": {
            "name": ["N", "NE", "NH1", "NH2"],
            "trail": ["CA", "CZ", "CZ", "CZ"],
            "element": ["N", "N", "N", "N"],
            "rotatable": [False, False, False, False],
            "state": ["N", "N", "N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        }
    },
    "LYS": {
        "donors": {
            "name": ["N", "NZ"],
            "trail": ["CA", "CE"],
            "element": ["N", "N"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        }
    },
    "HIS": {
        "donors": {
            "name": ["N", "ND1", "NE2"],
            "trail": ["CA", "CG", "CD2"],
            "element": ["N", "N", "N"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "N"]
        },
        "acceptors": {
            "name": ["O", "ND1", "NE2"],
            "element": ["O", "N", "N"],
            "state": ["N", "N", "N"]
        }
    },
    "TRP": {
        "donors": {
            "name": ["N", "NE1"],
            "trail": ["CA", "CD1"],
            "element": ["N", "N"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        }
    },
    "ASP": {
        "donors": {
            "name": ["N"],
            "trail": ["CA"],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "OD1", "OD2"],
            "element": ["O", "O", "O"],
            "state": ["N", "N", "N"]
        }
    },
    "GLU": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "OE1", "OE2"],
            "element": ["O", "O", "O"],
            "state": ["N", "N", "N"]
        }
    },
    "SER": {
        "donors": {
            "name": ["N", "OG"],
            "trail": [],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "O"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OG"],
            "element": ["O", "O"],
            "state": ["N", "N"]
        },
    },
    "CYS": {
        "donors": {
            "name": ["N", "SG"],
            "trail": [],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "S"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "SG"],
            "element": ["O", "S"],
            "state": ["N", "N"]
        },
    },
    "MET": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O", "SD"],
            "element": ["O", "S"],
            "state": ["N", "N"]
        },
    },
    "THR": {
        "donors": {
            "name": ["N", "OG1"],
            "trail": [],
            "trail": [["CA", "CB"], ["CB", "CA"]],
            "element": ["N", "O"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OG1"],
            "element": ["O", "O"],
            "state": ["N", "N"]
        },
    },
    "ASN": {
        "donors": {
            "name": ["N", "ND2"],
            "trail": [],
            "trail": [["CA", "CB"], ["CG", "CB"]],
            "element": ["N", "N"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OD1"],
            "element": ["O", "O"],
            "state": ["N", "N"]
        },
    },
    "GLN": {
        "donors": {
            "name": ["N", "NE2"],
            "trail": [],
            "trail": [["CA", "CB"], ["CD", "CG"]],
            "element": ["N", "N"],
            "rotatable": [False, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OE1"],
            "element": ["O", "O"],
            "state": ["N", "N"]
        },
    },
    "ALA": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "GLY": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "C"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "ILE": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "LEU": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "PRO": {
        "donors": {
            "name": [],
            "trail": [],
            "trail": [],
            "element": [],
            "rotatable": [],
            "state": []
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "PHE": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "TYR": {
        "donors": {
            "name": ["N", "OH"],
            "trail": [],
            "trail": [["CA", "CB"], ["CZ", "CE1"]],
            "element": ["N", "O"],
            "rotatable": [False, True],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["O", "OH"],
            "element": ["O", "O"],
            "state": ["N", "N"]
        },
    },
    "VAL": {
        "donors": {
            "name": ["N"],
            "trail": [],
            "trail": [["CA", "CB"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["O"],
            "element": ["O"],
            "state": ["N"]
        },
    },
    "G": {
        "donors": {
            "name": ["O2'", "N1", "N2", "O6"],
            "trail": [],
            "trail": [["C2'", "C1'"], ["C2", "N3"], ["C2", "N3"], ["C6", "C5"]],
            "element": ["O", "N", "N", "O"],
            "rotatable": [True, False, False, False],
            "state": ["N", "N", "N", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N3", "O6", "N7", "N1"],
            "element": ["O", "O", "O", "O", "O", "O", "N", "O", "N", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "A": {
        "donors": {
            "name": ["O2'", "N6"],
            "trail": [],
            "trail": [["C2'", "C1'"], ["C6", "C5"]],
            "element": ["O", "N"],
            "rotatable": [True, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N1", "N3", "N7"],
            "element": ["O", "O", "O", "O", "O", "O", "N", "N", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "C": {
        "donors": {
            "name": ["O2'", "N4"],
            "trail": [],
            "trail": [["C2'", "C1'"], ["C4", "C5"]],
            "element": ["O", "N"],
            "rotatable": [True, False],
            "state": ["N", "N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "N3", "O2"],
            "element": ["O", "O", "O", "O", "O", "O", "N", "O"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "U": {
        "donors": {
            "name": ["O2'", "N3", "O2", "O4"],
            "trail": [],
            "trail": [["C2'", "C1'"], ["C2", "N1"], ["C2", "N1"], ["C4", "N3"]],
            "element": ["O", "N", "O", "O"],
            "rotatable": [True, False, False, False],
            "state": ["N", "N", "T", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O2'", "O4", "O2", "N3"],
            "element": ["O", "O", "O", "O", "O", "O", "O", "O", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "DG": {
        "donors": {
            "name": ["N1", "N2", "O6"],
            "trail": [],
            "trail": [["C2", "N3"], ["C2", "N3"], ["C6", "C5"]],
            "element": ["N", "N", "O"],
            "rotatable": [False, False, False],
            "state": ["N", "N", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N3", "O6", "N7", "N1"],
            "element": ["O", "O", "O", "O", "O", "N", "O", "N", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N", "T"]
        }
    },
    "DA": {
        "donors": {
            "name": ["N6"],
            "trail": [],
            "trail": [["C6", "C5"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N1", "N3", "N7"],
            "element": ["O", "O", "O", "O", "O", "N", "N", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "DC": {
        "donors": {
            "name": ["N4"],
            "trail": [],
            "trail": [["C4", "C5"]],
            "element": ["N"],
            "rotatable": [False],
            "state": ["N"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "N3", "O2"],
            "element": ["O", "O", "O", "O", "O", "N", "O"],
            "state": ["N", "N", "N", "N", "N", "N", "N"]
        }
    },
    "DT": {
        "donors": {
            "name": ["N3", "O2", "O4"],
            "trail": [],
            "trail": [["C2", "N1"], ["C2", "N1"], ["C4", "N3"]],
            "element": ["N", "O", "O"],
            "rotatable": [False, False, False],
            "state": ["N", "T", "T"]
        },
        "acceptors": {
            "name": ["OP1", "OP2", "O5'", "O4'", "O3'", "O4", "O2", "N3"],
            "element": ["O", "O", "O", "O", "O", "O", "O", "N"],
            "state": ["N", "N", "N", "N", "N", "N", "N", "T"]
        }
    }
}
