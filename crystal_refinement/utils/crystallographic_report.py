# read cif and populate crystallographic data
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.text import WD_PARAGRAPH_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import matplotlib.pyplot as plt
from math import ceil
import numpy as np
import pandas as pd
from docx.shared import RGBColor
from docx.shared import Inches

import math
import argparse
from crystal_refinement.utils.cif_helper import *
from crystal_refinement.utils.crystallography_helper import *
from crystal_refinement.utils.space_groups import *



def black(heading):
    for run in heading.runs:
        run.font.color.rgb = RGBColor(0, 0, 0)


def add_mixed_text(paragraph, parts, clear=True):
    if clear: paragraph.clear()
    # print(parts)
    for part, formatting in parts:

        if 'bar' in formatting:
            bar_run = paragraph.add_run(f'{part}')

            bar_run = paragraph.add_run("\u0305")
            bar_run.font.superscript = True
            return
        run = paragraph.add_run(part)

        if 'it'in formatting:
            run.italic = True
        if 'sub' in formatting:
            run.font.subscript = True
        if 'bold' in formatting:
            run.font.bold = True
        if 'sup' in formatting:
            run.font.superscript = True
        if 'grey' in formatting:
            run.font.color.rgb = RGBColor(128, 128, 128)


def format_scientific(value, uncertainty, precision=4):


    if value == '':
        return ''
    
    # Handle uncertainty objects (from the uncertainties package)
    if hasattr(value, 'n'):
        val_num = value.n
        unc_num = value.std_dev
    else:
        val_num = float(value)
        unc_num = float(uncertainty)

    if math.isnan(unc_num) or unc_num == 0:
        return f"{val_num:g}"

    # Calculate significant figures for the uncertainty
    unc_scale = 10**np.floor(np.log10(unc_num))
    unc_digits = round(unc_num / unc_scale)
    
    # If rounding makes the uncertainty 10, adjust scale
    if unc_digits == 10:
        unc_scale *= 10
        unc_digits = 1

    formatted_val = f"{val_num:.{precision}f}"

    return f"{formatted_val}({int(unc_digits)})"


def format_coordinate(val, sigma, precision=4):

    if sigma == 0:
        if np.isclose(val, 0.5, atol=1e-4): return "1/2"
        if np.isclose(val, 0.25, atol=1e-4): return "1/4"
        if np.isclose(val, 0.75, atol=1e-4): return "3/4"
        if np.isclose(val, 0.2, atol=1e-4): return "1/5"
        if np.isclose(val, 0.4, atol=1e-4): return "2/5"
        if np.isclose(val, 0.6, atol=1e-4): return "3/5"
        if np.isclose(val, 0.8, atol=1e-4): return "4/5"
        if np.isclose(val, 0.0, atol=1e-4): return "0"
        if np.isclose(val, 1.0, atol=1e-4): return "1"
        if np.isclose(val, 0.3333, atol=1e-4): return "1/3"
        if np.isclose(val, 0.6667, atol=1e-4): return "2/3"

        return f"{round(val, precision)}"
    
    return format_scientific(val, sigma, precision=precision)


def create_word_table(refinement_data):

    space_group_nums_params = {}

    # triclinic
    for i in range(1, 3):
        space_group_nums_params[i] = "A-B-C_a-b-g"

    # monoclinic
    for i in range(3, 16):
        space_group_nums_params[i] = "A-B-C_b"

    # orthorhombic
    for i in range(16, 75):
        space_group_nums_params[i] = "A-B-C"

    # tetragonal
    for i in range(75, 143):
        space_group_nums_params[i] = "A-C"

    # trigonal
    for i in range(143, 168):
        space_group_nums_params[i] = "A-C"

    # hexagonal
    for i in range(168, 195):
        space_group_nums_params[i] = "A-C"

    # cubic
    for i in range(195, 231):
        space_group_nums_params[i] = "A"

    params = []
    sg_params = space_group_nums_params[int(refinement_data['_space_group_IT_number'])]

    if "_" in sg_params:
        params, angles = sg_params.split("_")
        params = params.split("-") if "-" in params else [params]
        angles = angles.split("-") if "-" in angles else [angles]
    else:
        params = sg_params.split("-") if "-" in sg_params else [sg_params]
        angles = []

    param_map = {
        "A": "a (Å)",
        "B": "b (Å)",
        "C": "c (Å)",
        "a": "alpha",
        "b": "beta",
        "g": "gamma",
    }

    # print(refinement_data.keys())
    # Data dictionary
    data = {
        "Formula": refinement_data["_chemical_formula_sum"].replace(" ", ""),
        "Formula mass (amu)": refinement_data["_chemical_formula_weight"],
        "Space group": refinement_data['_space_group_IT_number'].replace(" ", ""),
        "a (Å)": f"{refinement_data['_cell_length_a']}",
        "b (Å)": refinement_data['_cell_length_b'],
        "c (Å)": refinement_data['_cell_length_c'],
        "alpha": refinement_data['_cell_angle_alpha'],
        "beta": refinement_data['_cell_angle_beta'],
        "gamma": refinement_data['_cell_angle_gamma'],
        "V (Å³)": refinement_data['_cell_volume'],
        "Z": refinement_data['_cell_formula_units_Z'],
        "T (K)": refinement_data['_diffrn_ambient_temperature'],
        "Calculated density (g cm⁻³)": refinement_data['_exptl_crystal_density_diffrn'],
        
        "Crystal dimensions (mm)": f"{refinement_data['_exptl_crystal_size_max']} × {refinement_data['_exptl_crystal_size_mid']} × {refinement_data['_exptl_crystal_size_min']}",
        
        f"μ ({refinement_data['_diffrn_radiation_type']}) (mm⁻¹)".replace("\\", ""): refinement_data['_exptl_absorpt_coefficient_mu'],
        "2θ limits": f"{float(refinement_data.get('_diffrn_reflns_theta_min', 0.0))*2:.2f}–{float(refinement_data.get('_diffrn_reflns_theta_max', 0.0))*2:.2f}",
        "No. of data collected": refinement_data['_reflns_number_total'],
       
        "No. of unique data (Fo² < 0 included)": refinement_data['_reflns_number_total'],
        "No. of observed data (Fo² > 2σ(Fo²))": refinement_data['_reflns_number_gt'],
       
        "No. of variables": refinement_data['_refine_ls_number_parameters'],
       
        "Final R1a indices (all data)": f"{refinement_data['_refine_ls_R_factor_gt']} ({refinement_data['_refine_ls_R_factor_all']})",
        "Weighted wR2b factor (all data)": f"{refinement_data['_refine_ls_wR_factor_gt']} ({refinement_data['_refine_ls_wR_factor_ref']})",
        "Goodness of fit": f"{refinement_data['_refine_ls_goodness_of_fit_ref']}",
        "(Δρ)max, (Δρ)min (e Å⁻³)": f"{refinement_data['_refine_diff_density_max']}, {refinement_data['_refine_diff_density_min']}"
    }

    sample_name = data.get('Formula', 'Unknown Sample')

    # Create document
    doc = Document()
    
    sections = doc.sections
    if len(sections) > 0:
        section = sections[0]
        section.left_margin = Inches(1.4)

    doc.add_heading("", level=1)
    doc.add_paragraph(f"Table #. Crystallographic data from {sample_name}")

    for k in set(param_map.keys()) - set(params):
        data.pop(param_map[k])

    if "?" in data["Crystal dimensions (mm)"]:
        data.pop("Crystal dimensions (mm)")


    # Create table
    table = doc.add_table(rows=len(data), cols=2)
    table.autofit = True
    table.style = "Table Grid"

    italics_dict = {
        'Formula': [['Formula', 'bold']],
        'a (Å)': [['a', 'it'],[' (Å)', '']], 
        'b (Å)': [['b', 'it'],[' (Å)', '']], 
        'c (Å)': [['c', 'it'],[' (Å)', '']],
        'V (Å³)': [['V', 'it'], [' (Å', ''], ['3', 'sup'], [')', '']],
        'Z': [['Z', 'it']],
        'T (K)': [['T', 'it'], [' (K)', '']],
        'Calculated density (g cm⁻³)': [['Calculated density (g cm', ''], ['-3', 'sup'], [')', '']],
        'μ (MoKa) (mm⁻¹)': [['μ', 'it'], [' (Mo ', ''], ['K', 'it'], ['α', 'sub'], [ ') (mm', ''], ['-1', 'sup'], [')', '']],
        '2θ limits': [['2', ''], ['θ', 'it'], [' limits', ''], [' (\u00b0)', '']],
        'No. of unique data (Fo² < 0 included)': [['No. of unique data (', ''], ['F', 'it'], ['o', 'sub'], ['2', 'sup'], [' < 0 included)', '']],
        'No. of observed data (Fo² > 2σ(Fo²))': [['No. of observed data (', ''], ['F', 'it'], ['o', 'sub'], ['2', 'sup'], [' > 2', ''], ['σ', 'it'], ['(', ''], ['F', 'it'], ['o', 'sub'], ['2', 'sup'], ['))', '']],
        'Final R1a indices (all data)': [['Final ', ''], ['R', 'it'], ['1', 'sub'], [' ', ''], ['a', 'sup-it'], [' indices (all data)', '']],
        'Weighted wR2b factor (all data)': [['Weighted ', ''], ['wR', 'it'], ['2', 'sub'], ['b', 'sup-it'], [' factor (all data)', '']],
        '(Δρ)max, (Δρ)min (e Å⁻³)': [['(Δ', ''], ['ρ', 'it'], [')max, (Δ', ''], ['ρ', 'it'], [')min (e Å', ''], ['-3', 'sup'], [')', '']]
    }

    for i, (key, value) in enumerate(data.items()):
        row = table.rows[i]
        # Apply italics to keys as requested (a, c, V, Z, T, mu, etc)
        if key in italics_dict:
            cell_para = row.cells[0].paragraphs[0]
            add_mixed_text(cell_para, italics_dict.get(key))
        else:
            row.cells[0].paragraphs[0].add_run(key)

        # if f"{key}_us" in refinement_data:
        #     print(f"{key}_us", refinement_data[f"{key}_us"])
        
        if key == "Formula":
            cell_para = row.cells[1].paragraphs[0]
            add_mixed_text(cell_para, [[str(value), 'bold']])
        elif key == "Space group":
            sg = int(value)
            sg_latex = space_groups[sg]
            omml_element = latex_to_omml(sg_latex)

            para = row.cells[1].paragraphs[0]._element
            para.append(omml_element)
        else:
            cell_para = row.cells[1].paragraphs[0]
            cell_para.text = str(value)
        cell_para.alignment = WD_ALIGN_PARAGRAPH.CENTER

    p = doc.add_paragraph()
    add_mixed_text(p, [['a ', 'sup'], ['R', 'it'], ['1', 'sub'], ['(', ''], ['F', 'it'], [') = Σ||', ''], ['F', 'it'], ['o', 'sub'], ['| – |', ''], ['F', 'it'], ['c', 'sub'], ['|| / Σ|', ''], ['F', 'it'], ['o', 'sub'], ['|', '']])
    # p = doc.add_paragraph()
    add_mixed_text(p, [('\n', ''), ['b ', 'sup'], ['wR', 'it'], ['2', 'sup'], ['(', ''], ['F', 'it'], [') = [', ''], ['S', 'it'], ['[', ''], ['w', 'it'], ['(', ''], ['F', 'it'], ['o', 'sub'], ['2', 'sup'], [' - ', ''], ['F', 'it'], ['c', 'sub'], ['2', 'sup'], [')', ''], ['2', 'sup'], ['/', ''], ['S', 'it'], ['[', ''], ['w', 'it'], ['(', ''], ('F', 'it'), ('o', 'sub'), ('2', 'sup'), (')', ''), ('2', 'sup'), (']]', ''), ('1/2', 'sup'), 
                       
                       (' [', ''), ('w', 'it'), (' -1', 'sup'), (' = ', ''), ('σ', 'it'), ('2', 'sup'), ('(', ''), ('F', 'it'), ('o', 'sup'), (')', ''), ('2', 'sup'), (' + (0.0534', ''), ('P', 'it'), (')', ''), ('2', 'sup'), ('], where ', ''), ('P', 'it'), (' = ', ''), ('(', ''), ('F', 'it'), ('o', 'sub'), ('2', 'sup'), ('+2', ''), ('F', 'it'), ('c', 'sub'), ('2', 'sup'), (')/3', '')], clear=False)

    # Site table
    doc.add_paragraph(f"\nTable #. Atomic Coordinates and Equivalent Isotropic Displacement Parameters of {sample_name}")
    
    # sites, _ = get_wyckoff_symbol(refinement_data)
    sites = add_cifkit_labels(refinement_data)
    

    columns = [[('site', 'bold')],[('Wyckoff position', 'bold')], [('x', 'it-bold')], [('y', 'it-bold')], [('z', 'it-bold')], [('occupancy', 'bold')], [('U', 'bold-it'), ('eq', 'sub-bold'), ('a', 'sup'), (' (Å', 'bold'), ('2', 'sup-bold'), (')', 'bold')]]
    table = doc.add_table(rows=len(sites) + 1, cols=len(columns))
    table.style = 'Table Grid'

    header_row = table.rows[0]
    for col_idx, header_text in enumerate(columns):
        add_mixed_text(header_row.cells[col_idx].paragraphs[0], header_text)

    
    for i, site in enumerate(sites):
        row = table.rows[i + 1]
        row.cells[0].text = site.get('label', '')
        add_mixed_text(row.cells[1].paragraphs[0], [(str(site.get('multiplicity', '')), ''), (site.get('Wyckoff_symbol', ''), 'it')])
        # row.cells[2].text = f"{format_coordinate(*site.get('x', ''), precision=4)}"
        # row.cells[3].text = f"{format_coordinate(*site.get('y', ''), precision=4)}"
        # row.cells[4].text = f"{format_coordinate(*site.get('z', ''), precision=4)}"
        row.cells[2].text = f"{site.get('x', '').tv}"
        row.cells[3].text = f"{site.get('y', '').tv}"
        row.cells[4].text = f"{site.get('z', '').tv}"
        row.cells[5].text = f"{site.get('occupancy', '')}"
        row.cells[6].text = f"{site.get('Uiso', '')}"

        # z = site.get('x', '')
        # print(z.raw, z.nx, z.tv, z.n, z.num_dec, z.u, z.num_us)

    p = doc.add_paragraph()
    add_mixed_text(p, [('a', 'sup'), ('U', 'it'), ('eq', 'sub'), (' is defined as one-third of the trace of the orthogonalized ', ''), ('U', 'it'), ('ij', 'sub-it'), (' tensor.', '')])

    # coordination table
    cell_param_keys = [
        "_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"
    ]

    cell_params = [refinement_data[k] for k in cell_param_keys]
    # sigma_cell_params = [refinement_data[f"{k}_us"] for k in cell_param_keys]
    coord_table = get_coordination_data(refinement_data["cif_path"], sites, cell_params)

    
    coord_table_for_writing = []
    for ce in coord_table:
        center = ce[0]
        center_cn = ce[1]
        csum = 0

        for i, ne in enumerate(ce[2]):
            row = {'atom': '', 'neigh': ne[0], 'd': ne[1], 'c': ne[2], 'delta': ne[3], 'it': 0}
            csum += sum(row['c']) if isinstance(row['c'], list) else row['c']

            if i == 0:
                row['atom'] = center
            elif i == 1:
                row['atom'] = f"CN {center_cn}"

            if csum > center_cn:
                row['it'] = 1

            coord_table_for_writing.append(row)
        
        row = {'atom': '', 'neigh': '', 'd': '', 'c': '', 'delta': '', 'it': 0}
        coord_table_for_writing.append(row)

    if len(coord_table_for_writing) % 2 == 1:
        coord_table_for_writing = coord_table_for_writing[:-1]
    total_entries = len(coord_table_for_writing)
    # split_point = ceil(total_entries / 2)

    # Coord table
    # doc.add_paragraph(f"\nTable #. Coordination Environments in {sample_name}")
    p = doc.add_paragraph()
    add_mixed_text(p, [('Table #.Interatomic distances (', ''), ('d', 'it'), (', Å', ''), (')', ''), (', Δ values (Δ = 100(', ''), ('d', 'it'), ('Σ-', ''), ('r', 'it'), (')', ''), ('/Σ', ''), ('r', 'it'), (' is the sum of the respective atomic radii', ''), (' and atomic coordination numbers (CN) for ', ''), (f'{sample_name}.', '')])
    table = doc.add_table(rows=total_entries + 1, cols=5)
    table.style = 'Table Grid'
    table.autofit = True
    
    header_row = table.rows[0]
    header_row.cells[0].text = f"atom"
    add_mixed_text(header_row.cells[0].paragraphs[0], [['atom', 'bold']])
    header_row.cells[1].text = " "
    add_mixed_text(header_row.cells[2].paragraphs[0], [['d', 'it-bold'], [' (Å)', 'bold']])
    add_mixed_text(header_row.cells[3].paragraphs[0], [['count', 'bold']])
    add_mixed_text(header_row.cells[4].paragraphs[0], [['Δ', 'it-bold'], [' (%)', 'bold']])

    for i in range(5):
        header_row.cells[i].bold = True


    for i, entry in enumerate(coord_table_for_writing):
        target_row_idx = i + 1
        col_offset = 0 # Left side
        row = table.rows[target_row_idx]
        row.cells[0 + col_offset].text = str(entry['atom'])
        row.cells[1 + col_offset].text = str(entry['neigh'])
        if isinstance(entry['d'], list):
            # val = '/'.join([format_scientific(v, 0, precision=4) for v in entry['d']])
            val = format_scientific(entry['d'][0], 0, precision=4)
        else:
            val = format_scientific(entry['d'], 0, precision=4)
        row.cells[2 + col_offset].text = val
        if entry['c'] != '':
            val = int(sum(entry['c'])/len(entry['neigh'].split('/'))) if isinstance(entry['c'], list) else entry['c']
            row.cells[3 + col_offset].text = f"{val}\u00d7"
        
        if isinstance(entry['delta'], list):
            val = entry['delta']
            val = '/'.join([format_scientific(v, 0, precision=2) for v in val])
        else:
            val = format_scientific(entry['delta'], 0, precision=2)
        row.cells[4 + col_offset].text = val

        if entry.get('it', False):
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for run in paragraph.runs:
                        # run.italic = True
                        run.font.color.rgb = RGBColor(128, 128, 128)
    
    p = doc.add_paragraph()
    add_mixed_text(p, [('Neighbors outside of the coordination number are ', ''), ('in grey', 'grey'), ('.', '')])
    
    for para in doc.paragraphs:
        para.paragraph_format.first_line_indent = 0

    # Save file
    doc.save(f"crystal_structure_data_{sample_name}.docx")

    print("Word file created: crystal_structure_data.docx")


def cli_create_word_table():
    """CLI entry point for creating Word tables from CIF files."""
    parser = argparse.ArgumentParser(description="Generate a crystallographic report in Word format.")
    parser.add_argument("cif_path", help="Path to the CIF file")
    
    args = parser.parse_args()
    
    print(f"Processing CIF: {args.cif_path}")
    result = extract_cif_data(args.cif_path)
    create_word_table(result)


if __name__ == "__main__":
    cif_path = "/home/bala/Documents/44_shelxl/test_files/4_DyIrSn_done/test.cif"
    # cif_path = "/home/bala/Documents/44_shelxl/test_files/6_done_cell_moved/test.cif"
    # cif_path = "1534822.cif"
    
    # create_word_table(result)
    # cif_path = "/home/bala/Documents/data/not_prototype_CIFs/1232634.cif"
    result = extract_cif_data(cif_path)


    create_word_table(result)