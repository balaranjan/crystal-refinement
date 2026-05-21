from crystal_refinement.utils.crystallography_helper import CFloat
from crystal_refinement.SHELX import SHELXFile
from fractions import Fraction
import numpy as np
import subprocess
import argparse
import sys
import re
import os
import traceback


def run_platon(platon_executable, cif_path):

    if not os.path.isfile(cif_path):
        raise FileNotFoundError(f"Input file not found: {cif_path}")

    try:
        # check=True ensures that a CalledProcessError is raised if the exit code is non-zero
        process = subprocess.Popen(
            [platon_executable, cif_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1           # Decode output to string
        )
        process.stdin.write("stidy" + '\n')
        process.stdin.flush() 
        # time.sleep(0.1)

        process.stdin.write("quit" + '\n')
        process.stdin.flush()
        # time.sleep(0.1)

        # process.stdin.close()
        stdout, stderr = process.communicate(timeout=60)

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, platon_executable, output=stdout, stderr=stderr)

        return stdout

    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(f"Error Output: {e.stderr}")
        raise
    except Exception as e:
        print(f"System error occurred: {e}")
        raise


def run_stidy(platon_executable, cif_path):

    sty_path = f"{cif_path[:-3]}sty"

    if not os.path.isfile(sty_path):
        raise FileNotFoundError(f"Input file not found: {sty_path}")
    try:
        result = subprocess.run(
            [platon_executable, '-Y', sty_path],
            check=True,           
            capture_output=True,  
            text=True              
        )

        return result.stdout

    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print(f"Error Output: {e.stderr}")
        raise
    except Exception as e:
        print(f"System error occurred: {e}")
        raise


def extract_platon_data(text):
    """
    Extracts unit cell parameters and site data from Structure Tidy output.
    """
    data = {
        'a': 0.0, 'b': 0.0, 'c': 0.0, 
        'alpha': 0.0, 'beta': 0.0, 'gamma': 0.0,
        'site_data': []
    }

    cell_pattern = r"Cell\s*:\s*([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)"
    cell_match = re.search(cell_pattern, text)
    
    if cell_match:
        data['a'] = CFloat(cell_match.group(1))
        data['b'] = CFloat(cell_match.group(2))
        data['c'] = CFloat(cell_match.group(3))
        data['alpha'] = CFloat(cell_match.group(4))
        data['beta'] = CFloat(cell_match.group(5))
        data['gamma'] = CFloat(cell_match.group(6))

    # Extract Site Data
    lines = text.splitlines()
    
    def parse_site(line):
        line = line.split()
        if len(line) < 6 or line[0] in ['Wyckoff', 'Volume']:
            return 
        if not line[-1].isnumeric():
            line.append('')

        label = f"{line[-2]}{line[-1]}"
        vals = [label, line[1], *[CFloat(v, no_frac=True) for v in line[2:5]]]
        if len(line) == 7:
            vals.append(1)
            vals.extend(line[5:])
        else:
            vals.append(CFloat(line[5]))
            vals.extend(line[6:])
        return vals
    
    for line in lines[6:]:
        vals = parse_site(line)
        if vals:
            data['site_data'].append(vals)

    return data


def update_cif(cif_text, params):

    def parse_wyckoff(s):
        num = s[:s.index('(')]
        letter = s[s.index('(')+1:s.index(')')]
        return num, letter

    def get_unc(token):
        return token[token.index('('):] if '(' in token else ''

    lines = cif_text.splitlines()

    old_pos_unc = {}
    other_data = {}
    in_site_loop = False
    site_col_indices = {}
    for i in range(len(lines)):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith('loop_'):
            in_site_loop = False
            site_col_indices = {}
        elif in_site_loop and stripped.startswith('_atom_site'):
            tag = stripped.split()[0]
            site_col_indices[tag] = len(site_col_indices)
        elif in_site_loop and stripped and not stripped.startswith('_') and not stripped.startswith('loop_'):
            tokens = stripped.split()
            other_data[tokens[0]] = [CFloat(tokens[5]), CFloat(tokens[7])]
            if tokens and not tokens[0].startswith('_'):
                lbl = tokens[0]
                xi = site_col_indices.get('_atom_site_fract_x')
                yi = site_col_indices.get('_atom_site_fract_y')
                zi = site_col_indices.get('_atom_site_fract_z')
                old_pos_unc[lbl] = (
                    get_unc(tokens[xi]) if xi is not None and xi < len(tokens) else '',
                    get_unc(tokens[yi]) if yi is not None and yi < len(tokens) else '',
                    get_unc(tokens[zi]) if zi is not None and zi < len(tokens) else '',
                )
        if '_atom_site_label' in stripped and '_geom' not in stripped:
            in_site_loop = True
            site_col_indices = {'_atom_site_label': 0}

    result = []
    i = 0
    i_end = len(lines)

    while i < len(lines):
        line = lines[i]

        # Update cell lengths
        if line.strip().startswith('_cell_length_a'):
            old_val = line.split()[1]; 
            unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['a'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['a']))
        elif line.strip().startswith('_cell_length_b'):
            old_val = line.split()[1]; unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['b'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['b']))
        elif line.strip().startswith('_cell_length_c'):
            old_val = line.split()[1]; unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['c'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['c']))

        # Update cell angles
        elif line.strip().startswith('_cell_angle_alpha'):
            old_val = line.split()[1]; unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['alpha'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['alpha']))
        elif line.strip().startswith('_cell_angle_beta'):
            old_val = line.split()[1]; unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['beta'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['beta']))
        elif line.strip().startswith('_cell_angle_gamma'):
            old_val = line.split()[1]; unc = old_val[old_val.index('('):] if '(' in old_val else ''
            params['gamma'].update_u(unc)
            result.append(line[:line.index(line.split()[1])] + str(params['gamma']))

        # Replace atom_site loop
        elif line.strip().startswith('loop_') and i + 1 < len(lines) and '_atom_site_label' in lines[i+1] and 'geom' not in lines[i+1]:
            # Skip old loop header and data
            result.append('loop_')
            i += 1
            while i < len(lines) and (lines[i].strip().startswith('_atom_site') or lines[i].strip() == ''):
                i += 1

            # Write new loop headers
            result.append(' _atom_site_label')
            result.append(' _atom_site_type_symbol')
            result.append(' _atom_site_symmetry_multiplicity')
            result.append(' _atom_site_Wyckoff_symbol')
            result.append(' _atom_site_fract_x')
            result.append(' _atom_site_fract_y')
            result.append(' _atom_site_fract_z')
            result.append(' _atom_site_U_iso_or_equiv')
            result.append(' _atom_site_occupancy')

            # Write new atom site data
            for site in params['site_data']:
                # print(site)
                label, wyck_str, x, y, z, occ, sym, _ = site
                mult, letter = parse_wyckoff(wyck_str)
                ux, uy, uz = old_pos_unc.get(label, ('', '', ''))
                x.update_u(ux)
                y.update_u(uy)
                z.update_u(uz)

                result.append(f"{label} {sym} {mult} {letter} {x} {y} {z} {other_data[label][0]} {other_data[label][1]}")
                i += 1
            i_end = i

            # Skip old atom data lines
            while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('loop_') and not lines[i].strip().startswith('_'):
                i += 1
            continue
    
        else:
            result.append(line)
        
        i += 1

    output = result[:i_end]
    i_aniso = [i for i in range(len(output)) if "atom_site_aniso" in output[i]]

    if len(i_aniso):
        output = output[:i_aniso[0]-1]

    output.append("\n")
    for line in lines[i_end:]:
        if "_refine_diff_density" in line:
            output.append(line)

    return '\n'.join(output)


def standardize(cif_path):

    original_dir = os.getcwd()
    cif_path = cif_path.split(os.sep)
    cif_filename = cif_path[-1]
    cif_path = os.sep.join(cif_path[:-1])
    if cif_path == "":
        cif_path = "."

    try:
        os.chdir(cif_path)
        run_platon("platon", cif_filename)
        data = run_stidy('platon', cif_filename)
        # print(data)
        data = extract_platon_data(data)  

        updated_cif = update_cif(open(cif_filename, 'r').read(), data)
        with open(f"{cif_filename[:-4]}_standardized.cif", "w") as f:
            f.write(updated_cif)

    except:
        print("Error --")
        print(traceback.format_exc())
    finally:
        os.chdir(original_dir)


def cli_standardize():
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_path")
    parser.add_argument("--platon-executable", default="platon")
    
    args = parser.parse_args()
    
    try:
        standardize(args.cif_path)
        print("Standardization complete.")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        print(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    cif_path = "/home/bala/Documents/44_shelxl/test_files/6_done_cell_moved/test.cif"
    run_platon("platon", cif_path)
    data = run_stidy('platon', cif_path)
    data = extract_platon_data(data)

    for k, v in data.items():
        print(f"{k:<10} {v}")