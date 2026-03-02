from collections import defaultdict
import re


class Composition:
    def __init__(self, formula):
        self.formula = formula
        self.formula_dict = self._parse_formula(self.formula)

    def get_atomic_fraction(self, element):
        return self.formula_dict[element] / sum(self.formula_dict.values())
    
    def _parse_formula(self, formula: str, strict: bool = True) -> dict[str, float]:
        """Credits: pymatgen team.

        Args:
            formula (str): A string formula, e.g. Fe2O3, Li3Fe2(PO4)3.
            strict (bool): Whether to throw an error if formula string
            is invalid (e.g. empty).
                Defaults to True.

        Returns:
            Composition with that formula.

        Notes:
            In the case of Metallofullerene formula (e.g. Y3N@C80),
            the @ mark will be dropped and passed to parser.
        """
        # Raise error if formula contains special characters
        # or only spaces and/or numbers

        if "'" in formula:
            formula = formula.replace("'", "")

        if strict and re.match(r"[\s\d.*/]*$", formula):
            raise ValueError(f"Invalid formula={formula}")

        # For Metallofullerene like "Y3N@C80"
        formula = formula.replace("@", "")
        # Square brackets are used in formulas to denote coordination
        # complexes (gh-3583)

        formula = formula.replace("[", "(")
        formula = formula.replace("]", ")")

        def get_sym_dict(form: str, factor: float) -> dict[str, float]:
            sym_dict: dict[str, float] = defaultdict(float)
            for match in re.finditer(r"([A-Z][a-z]*)\s*([-*\.e\d]*)", form):
                el = match[1]
                amt = 1.0
                if match[2].strip() != "":
                    amt = float(match[2])
                sym_dict[el] += amt * factor
                form = form.replace(match.group(), "", 1)
            if form.strip():
                raise ValueError(f"{form} is an invalid formula!")
            return sym_dict

        match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        while match:
            factor = 1.0
            if match[2] != "":
                factor = float(match[2])
            unit_sym_dict = get_sym_dict(match[1], factor)
            expanded_sym = "".join(
                f"{el}{amt}" for el, amt in unit_sym_dict.items()
            )
            expanded_formula = formula.replace(match.group(), expanded_sym, 1)
            formula = expanded_formula
            match = re.search(r"\(([^\(\)]+)\)\s*([\.e\d]*)", formula)
        return get_sym_dict(formula, 1)


class Element:
    def __init__(self, name):
        if len(name) > 1:
            self.name = name[0]
            for c in name[1:]:
                self.name += c.lower()
            self.name = self.name.strip()
        else:
            self.name = name
        self.number = element_data[self.name][0]
        self.atomic_radius = element_data[self.name][1]


element_data = {
'H': [1, 0.25],
'He': [2, None],
'Li': [3, 1.45],
'Be': [4, 1.05],
'B': [5, 0.85],
'C': [6, 0.7],
'N': [7, 0.65],
'O': [8, 0.6],
'F': [9, 0.5],
'Ne': [10, None],
'Na': [11, 1.8],
'Mg': [12, 1.5],
'Al': [13, 1.25],
'Si': [14, 1.1],
'P': [15, 1.0],
'S': [16, 1.0],
'Cl': [17, 1.0],
'Ar': [18, 0.71],
'K': [19, 2.2],
'Ca': [20, 1.8],
'Sc': [21, 1.6],
'Ti': [22, 1.4],
'V': [23, 1.35],
'Cr': [24, 1.4],
'Mn': [25, 1.4],
'Fe': [26, 1.4],
'Co': [27, 1.35],
'Ni': [28, 1.35],
'Cu': [29, 1.35],
'Zn': [30, 1.35],
'Ga': [31, 1.3],
'Ge': [32, 1.25],
'As': [33, 1.15],
'Se': [34, 1.15],
'Br': [35, 1.15],
'Kr': [36, None],
'Rb': [37, 2.35],
'Sr': [38, 2.0],
'Y': [39, 1.8],
'Zr': [40, 1.55],
'Nb': [41, 1.45],
'Mo': [42, 1.45],
'Tc': [43, 1.35],
'Ru': [44, 1.3],
'Rh': [45, 1.35],
'Pd': [46, 1.4],
'Ag': [47, 1.6],
'Cd': [48, 1.55],
'In': [49, 1.55],
'Sn': [50, 1.45],
'Sb': [51, 1.45],
'Te': [52, 1.4],
'I': [53, 1.4],
'Xe': [54, None],
'Cs': [55, 2.6],
'Ba': [56, 2.15],
'La': [57, 1.95],
'Ce': [58, 1.85],
'Pr': [59, 1.85],
'Nd': [60, 1.85],
'Pm': [61, 1.85],
'Sm': [62, 1.85],
'Eu': [63, 1.85],
'Gd': [64, 1.8],
'Tb': [65, 1.75],
'Dy': [66, 1.75],
'Ho': [67, 1.75],
'Er': [68, 1.75],
'Tm': [69, 1.75],
'Yb': [70, 1.75],
'Lu': [71, 1.75],
'Hf': [72, 1.55],
'Ta': [73, 1.45],
'W': [74, 1.35],
'Re': [75, 1.35],
'Os': [76, 1.3],
'Ir': [77, 1.35],
'Pt': [78, 1.35],
'Au': [79, 1.35],
'Hg': [80, 1.5],
'Tl': [81, 1.9],
'Pb': [82, 1.8],
'Bi': [83, 1.6],
'Po': [84, 1.9],
'At': [85, None],
'Rn': [86, None],
'Fr': [87, None],
'Ra': [88, 2.15],
'Ac': [89, 1.95],
'Th': [90, 1.8],
'Pa': [91, 1.8],
'U': [92, 1.75],
'Np': [93, 1.75],
'Pu': [94, 1.75],
'Am': [95, 1.75],
'Cm': [96, None],
'Bk': [97, None],
'Cf': [98, None],
'Es': [99, None],
}