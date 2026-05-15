from crystal_refinement.utils.crystallography_helper import CFloat
from collections import defaultdict
import math
import re


def get_float(s):

    if "(" in s:
        s, us = s.split("(")
        num_zeros = len(s.split('.')[1])
        us = float(us.replace(")", "")) / math.pow(10, num_zeros)
        return float(s), us
    else:
        return float(s), 0.0
    

def extract_cif_data(filepath):

    data = {"cif_path": filepath}

    with open(filepath, 'r') as f:
        content = f.read()
        lines = content.splitlines()
        data['site_data'] = get_site_data(lines)
    
    # Helper: extract simple key-value pairs
    def get_value(key, text):
       
        pattern = rf'^{re.escape(key)}\s+(.+)$'
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            return match.group(1).strip().strip("'")
        # Value might be on the next line (e.g., _chemical_formula_sum)
        pattern2 = rf'^{re.escape(key)}\s*\n\s*\'?([^\n\']+)\'?'
        match2 = re.search(pattern2, text, re.MULTILINE)
        if match2:
            return match2.group(1).strip().strip("'")
        return None
    
    numeric = [
        '_cell_length_a',
        '_cell_length_b',
        '_cell_length_c',
        '_cell_angle_alpha',
        '_cell_angle_beta',
        '_cell_angle_gamma',
        '_cell_volume',
    ]

    cif_keys = [
        '_chemical_formula_sum',
        '_chemical_formula_weight',
        '_space_group_crystal_system',
        '_space_group_IT_number',
        '_space_group_name_H-M_alt',
        '_cell_formula_units_Z',
        '_cell_length_a',
        '_cell_length_b',
        '_cell_length_c',
        '_cell_angle_alpha',
        '_cell_angle_beta',
        '_cell_angle_gamma',
        '_cell_volume',
        '_exptl_crystal_density_diffrn',
        '_exptl_crystal_size_max',       
        '_exptl_crystal_size_mid',       
        '_exptl_crystal_size_min',       
        '_exptl_absorpt_coefficient_mu',
        '_exptl_crystal_F_000',
        '_diffrn_radiation_type',
        '_diffrn_ambient_temperature',
        '_diffrn_reflns_theta_min',
        '_diffrn_reflns_theta_max',
        '_diffrn_reflns_number',
        '_reflns_number_total',
        '_reflns_number_gt',
        '_refine_ls_number_parameters',    
        '_refine_ls_number_restraints',   
        '_refine_ls_R_factor_all',        
        '_refine_ls_R_factor_gt',         
        '_refine_ls_wR_factor_ref',       
        '_refine_ls_wR_factor_gt',        
        '_refine_ls_goodness_of_fit_ref', 
        '_refine_ls_restrained_S_all',    
        '_refine_diff_density_max',
        '_refine_diff_density_min'
    ]

    for key in cif_keys:
        value = get_value(key, content)
        if key in numeric:
            data[key] = CFloat(value)
            # data[f"{key}_str"] = value
            # value, us = get_float(value)
            # data[f"{key}_us"] = us

        data[key] = value

    return data



def _parse_formula(formula: str, strict: bool = True) -> dict[str, float]:
    """Copied from pymatgen.

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


def get_site_data(lines):

    headers_numeric = ["x", "y", "z", "occupancy", "multiplicity"]
    i_start, i_end = None, None

    for i, line in enumerate(lines):
        if "loop_" in line and "atom_site" in lines[i + 1]:
            i_start = i + 1
            break

    for i, line in enumerate(lines[i_start + 1 :]):
        if (
            line.startswith("_")
            or line.startswith("#")
            or line.startswith("loop_")
        ) and "atom_site" not in line:
            i_end = i + i_start + 1
            break

    values = lines[i_start:i_end]

    headers = [
        line.replace("_", " ")
        .strip()
        .replace("\n", "")
        .replace("atom site ", "")
        for line in values
        if "atom_site" in line
        and "atom_site_aniso" not in line
        and "atom_site_iso" not in line
    ]

    site_labels = defaultdict(int)
    site_data = []

    header_map = {
        "label": "label",
        "type symbol": "symbol",
        "symmetry multiplicity": "multiplicity",
        "Wyckoff symbol": "Wyckoff_symbol",
        "fract x": "x",
        "fract y": "y",
        "fract z": "z",
        "occupancy": "occupancy",
        "U iso or equiv": "Uiso"
    }

    # sig_figs = defaultdict(int)

    for site in values[len(headers) :]:
        site = site.replace("\n", "").strip().split()
        if not len(site):
            continue
        _site_values = {}
        for header, value in zip(headers, site):
            header = header_map.get(header)

            if not header:
                continue

            if header == "label":
                _site_values["label"] = value
            if header == "symbol":
                value = value.replace("-", "").replace("+", "")
                value = list(_parse_formula(value).keys())[0]
                site_labels[value] += 1

                _site_values[header] = value
            elif header in headers_numeric:
                # if header in ["x", "y", "z"] and "(" in value:
                #     sfs = value.split('.')[1].split('(')[0]
                    # sig_figs[len(sfs)] += 1
                
                value = CFloat(value)
                # if header == "multiplicity":
                #     value = int(value[0])
                _site_values[header] = value
            else:
                _site_values[header] = value

        if len(_site_values):
            site_data.append(_site_values)

    for i, v in enumerate(site_data):
        for k in ['x', 'y', 'z']:
            _val = v[k]
            _val.update_n(_val.n % 1.0)
        site_data[i] = v

    # if len(sig_figs):
    #     sig_figs = sorted([[k, v] for k, v in sig_figs.items()], key=lambda kv: kv[1])
    #     sig_figs = sig_figs[-1][0]
    # else:
    #     sig_figs = 4
    return site_data

