from crystal_refinement.utils.element_data import element_data
from collections import defaultdict
import re


class Composition:
    def __init__(self, formula):
        self.formula = formula
        self.formula_dict = dict(self._parse_formula(self.formula))
        self.formula_dict = dict(sorted(self.formula_dict.items(), key=lambda item: element_data[item[0]][2]))

    def get_atomic_fraction(self, element):
        if element in self.formula_dict:
            return self.formula_dict[element] / sum(self.formula_dict.values())
        else:
            return 0.0
    
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

