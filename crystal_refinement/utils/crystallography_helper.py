import gemmi
import spglib
import numpy as np
import pandas as pd
from crystal_refinement.utils.site_plots import CN_of_site, CN_of_site_stable
from crystal_refinement.utils.element_data import element_data
from uncertainties import ufloat, umath
from cifkit import Cif
from cifkit.utils import unit
from fractions import Fraction
from crystal_refinement.utils.cfloat import CFloat
import os



def fractional_to_cartesian(
    fractional_coords: list[float],
    cell_lengths: list[float],
    cell_angles_rad: list[float],
) -> list[float]:
    """Convert fractional coordinates to Cartesian coordinates using
    cell lengths and angles."""

    alpha, beta, gamma = cell_angles_rad
    alpha *= np.pi/180.0
    beta *= np.pi/180.0
    gamma *= np.pi/180.0

    # Calculate the components of the transformation matrix
    a, b, c = cell_lengths
    cos_alpha = umath.cos(alpha)
    cos_beta = umath.cos(beta)
    cos_gamma = umath.cos(gamma)
    sin_gamma = umath.sin(gamma)

    # The volume of the unit cell
    volume = (
        a
        * b
        * c
        * umath.sqrt(
            1
            - cos_alpha**2
            - cos_beta**2
            - cos_gamma**2
            + 2 * cos_alpha * cos_beta * cos_gamma
        )
    )

    # Transformation matrix from fractional to Cartesian coordinates
    matrix = np.array(
        [
            [a, b * cos_gamma, c * cos_beta],
            [
                0,
                b * sin_gamma,
                c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma,
            ],
            [0, 0, volume / (a * b * sin_gamma)],
        ]
    )

    cartesian_coords = np.dot(matrix, fractional_coords).flatten()

    return cartesian_coords



def add_cifkit_labels(data):
    cif = Cif(data['cif_path'])
    unitcell_points = cif.unitcell_points
    loop_vals = cif._loop_values
    site_symbol_map = dict(zip([l for l in loop_vals[0]], [s for s in loop_vals[1]]))

    site_data = data['site_data']
    site_data = {s['label']: s for s in site_data}

    # match
    for clabel in site_symbol_map.keys():
        for k, v in site_data.items():
            for x, y, z, cl in unitcell_points:
                
                if cl != clabel or site_symbol_map[cl] != v['symbol']:
                    continue

                x %= 1
                y %= 1
                z %= 1

                if np.allclose([v['x'].n, v['y'].n, v['z'].n], [x, y, z], atol=1e-4):
                    v['cifk_label'] = clabel
                    site_data[k] = v
                    break  

    return list(site_data.values())


def mixing_sites(unitcell_points):
    mixing = {}
    positions = np.array([p[:-1] for p in unitcell_points])
    labels = np.array([p[-1] for p in unitcell_points])

    for p, l in zip(positions, labels):
        d = np.linalg.norm(positions - p, axis=1).squeeze()
        same_pos = labels[d == 0]
        if len(same_pos) > 1:
            mixing[str(l)] = '/'.join(sorted(same_pos))

    return mixing




def get_coordination_data(cif_path, site_data, cell_params, nround=4):

    cif = Cif(cif_path)
    cif.compute_connections()
    cif.compute_CN()
    

    loop_vals = cif._loop_values
    site_symbol_map = dict(zip([l for l in loop_vals[0]], [s for s in loop_vals[1]]))

    conns = cif.connections
    table_data = []
    site_data = {s['cifk_label']: s for s in site_data}


    mixing_data = mixing_sites(cif.unitcell_points)

    _supercell_points = cif.supercell_points
    _supercell_points_u = []
    for x, y, z, label in _supercell_points:
        us = site_data[label]
        _x = us['x']; _x.update_n(x)
        _y = us['y']; _y.update_n(y)
        _z = us['z']; _z.update_n(z)
        _supercell_points_u.append([_x.uf, _y.uf, _z.uf, label])

    a, b, c, alpha, beta, gamma = [CFloat(v).uf for v in cell_params]
    lengths_u = [a, b, c]
    angles_u = [alpha, beta, gamma]

    lengths = cif.unitcell_lengths
    angles = cif.unitcell_angles

    supercell_points, supercell_points_u, supercell_labels = [], [], []
    for i in range(len(_supercell_points_u)):
        x, y, z, s = _supercell_points_u[i]
        cart_u = fractional_to_cartesian([x, y, z],
                                            lengths_u,
                                            angles_u)
        supercell_points_u.append(np.array(cart_u))

        x, y, z, s = _supercell_points[i]
        cart = unit.fractional_to_cartesian([x, y, z],
                                            lengths,
                                            angles)
        supercell_points.append(np.array(cart))
        supercell_labels.append(s)

    supercell_points_u = np.array(supercell_points_u)
    supercell_points = np.array(supercell_points)

    selected_CNs = None
    if os.path.isfile("CN.txt"):
        user_prefs = pd.read_csv("CN.txt")
        selected_CNs = dict(zip(user_prefs['Site'].tolist(), user_prefs['CN'].tolist()))
        print(selected_CNs)
        print(f"CN values from CN.txt will be used for sites lited in CN.txt")
    
    mixing_labels_calcd = []
    for site, neighbors in conns.items():
        if site in mixing_data:
            if mixing_data[site] in mixing_labels_calcd:
                continue
            mixing_labels_calcd.append(mixing_data[site])

        neighbors = sorted(neighbors, key=lambda x: x[1])[:21]
        
        if selected_CNs and site in selected_CNs:
            CN = selected_CNs[site]
        else:
            CN = CN_of_site(neighbors[:21], verbose=False)

        site_coord_w_u = supercell_points_u[np.argmin(np.linalg.norm(supercell_points - np.array(neighbors[0][2]), axis=1))]


        neighbors_table = []
        r1 = element_data[site_symbol_map[site]][1]
        for neighbor in neighbors[:20]:
            label = neighbor[0]
            neighbor_dists = np.linalg.norm(supercell_points - neighbor[3], axis=1).squeeze()
            ind_neighbor_coord_w_u = np.argmin(neighbor_dists)
            neighbor_coord_w_u = supercell_points_u[ind_neighbor_coord_w_u]

            

            dist = umath.sqrt(sum(umath.pow(site_coord_w_u[i] - neighbor_coord_w_u[i], 2)
             for i in range(3)))
            
            delta = "N.A."
            if label in mixing_data:
                labels = mixing_data[label].split('/')
                for _label in labels:
                    r2 = element_data[site_symbol_map[_label]][1]
                    if r1 and r2:
                        s = r1 + r2
                        delta = 100 * ((dist - s)/s)
                    neighbors_table.append({'neigh': _label, 'd': dist, 'delta': delta, 'd_num': round(dist.n, nround), 'count': 0, 'CN': CN})
            else:
                r2 = element_data[site_symbol_map[neighbor[0]]][1]
                if r1 and r2:
                    s = r1 + r2
                    delta = 100 * ((dist - s)/s)
                neighbors_table.append({'neigh': label, 'd': dist, 'delta': delta, 'd_num': round(dist.n, nround), 'count': 0, 'CN': CN})


        neighbors_table = neighbors_table
        neighbors_table = pd.DataFrame(neighbors_table)
        neighbors_table = (
            neighbors_table.groupby(['d_num', 'neigh'])
            .agg({
                # 'neigh': 'first',
                'd': 'first',
                'delta': 'first',
                'count': 'size',
                'CN': 'first'
            })
            .reset_index()
            .sort_values('d_num')
            .apply(lambda row: [row['neigh'], row['d'], row['count'], row['delta'], row['CN']], axis=1)
            .tolist()
        )

        if mixing_data:
            neighbors_table_mix = []
            for d in sorted(list(set([v[1].n for v in neighbors_table]))):
                same_d = [v for v in neighbors_table if v[1].n == d]
                if len(same_d) > 1:
                    same_inds = []
                    for i in range(len(same_d)):
                        if same_d[i][0] in mixing_data:
                            same_inds.append(i)
                            val = same_d[i]
                            val[0] = mixing_data[val[0]]
                            same_d[i] = val

                    merged_row = [same_d[same_inds[0]][0]]
                    keys = [same_d[i][0] for i in same_inds]
                    for j in [1, 2, 3]:
                        val = [same_d[k][j] for k in same_inds]
                        sd = dict(zip(val, keys))
                        val = sorted(val, key=lambda x: sd[x])
                        merged_row.append(val)
                    merged_row.append(same_d[same_inds[0]][-1])
                    neighbors_table_mix.append(merged_row)

                    neighbors_table_mix.extend([same_d[i] for i in range(len(same_d)) if i not in same_inds])
                else:
                    neighbors_table_mix.extend(same_d)
            neighbors_table = neighbors_table_mix

        if site in mixing_data:
            site = mixing_data[site]
        table_data.append([site, CN, neighbors_table])
        
    return table_data