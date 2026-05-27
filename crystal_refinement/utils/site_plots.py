import os
import argparse
from cifkit import Cif
from cifkit.utils import unit
import numpy as np
import pandas as pd
import pyvista as pv
from itertools import combinations
from scipy.spatial import ConvexHull
from crystal_refinement.utils.element_data import element_data
from crystal_refinement.utils.cfloat import CFloat
from matplotlib.colors import ListedColormap
import ast
pv.OFF_SCREEN = True
import traceback


def get_element_color(element):
    R = ["H", "Li", "Na", "K", "Rb", "Cs", "Fr",          # Group 1A
     "Be", "Mg", "Ca", "Sr", "Ba", "Ra",               # Group 2A
     "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",  # Lanthanides
     "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
     "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",  # Actinides
     "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Zr"]

    M = ["Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",  # Period 4 TM
        "Y",  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",  # Period 5 TM
        "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",         # Period 6 TM
        "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",]         # Period 7 TM

    X = ["B",  "Si", "Ge", "As", "Sb", "Te", "Po", "Sn", "Al"]                      # Metalloids

    if element in R:
        return "blue"
    elif element in M:
        return "grey"
    elif element in X:
        return "red"
    else:
        return "green"


def CN_of_site(v, full=False, verbose=False):
    
    # Finds the coordination numbers using the d/d_min method.
    points_wd =[[p[3], p[1], p[0]] for p in v]
    
    # sort
    points_wd = sorted(points_wd, key=lambda x: x[1])[:30]
    distances = np.array([p[1] for p in points_wd])
    distances /= distances.min()

    gaps = np.array([round(distances[i] - distances[i-1], 4) for i in range(1, len(distances))])
    ind_gaps = np.argsort(gaps, stable=True)

    CN_values = np.array(ind_gaps[::-1]) + 1
    CN_values = CN_values[CN_values >= 4]
    if verbose:
        print(gaps[ind_gaps[::-1]])
        print(CN_values)

    if full: return CN_values
    return CN_values[0]


def CN_of_site_stable(v, verbose=False):
    """
    Finds coordination numbers using a stabilized d/d_min method.
    Groups neighbors by distance tolerance to handle ties.
    """
    # Extract distances only (column index 1 in your tuple structure)
    distances = np.array([p[1] for p in v])
    
    if len(distances) == 0:
        return 0
    
    min_dist = distances.min()
    normalized_distances = distances / min_dist
    
    # Sort by normalized distance
    sorted_indices = np.argsort(normalized_distances)
    sorted_dists = normalized_distances[sorted_indices]
    
    gaps = np.diff(sorted_dists)
    
    # Handle ties: if gap is very small, it's not a real "jump"
    tolerance = 1e-4  
    significant_gaps = np.where(gaps > tolerance, gaps, 0)
    
    if verbose:
        print(f"Normalized distances: {sorted_dists}")
        print(f"Gaps: {gaps}")
        print(f"Significant gaps: {significant_gaps}")

    first_significant_gap_idx = np.argmax(significant_gaps)
    
    max_gap = significant_gaps.max()
    if max_gap < tolerance:
        CN = len(distances)
    else:
        CN = first_significant_gap_idx + 1
        
    if verbose:
        print(f"Coordination Number: {CN}")
        
    return CN


def get_points_and_labels(point, label, zmin=0.0, zmax=1.0):
    points = [point]
    labels = [label]

    if np.any(point == zmin):
        _point = point.copy()
        _point[_point==zmin] = zmax
        points.append(_point)
        labels.append(label)

    if np.any(point == zmax):
        _point = point.copy()
        _point[_point==zmin] = zmin
        points.append(_point)
        labels.append(label)

    return points, labels


def label_position(mesh, selected_axes, offset=0.5):
    """Returns a point slightly to the right of the mesh bounding box."""

    axis_vertical, axis_horizontal = selected_axes
    xmax = (mesh.bounds[axis_horizontal*2] + mesh.bounds[(axis_horizontal*2)+1]) / 2
    center_y = (mesh.bounds[axis_vertical*2] + mesh.bounds[(axis_vertical*2)+1]) / 2
    center_y = mesh.bounds[axis_vertical*2]

    iz = [i for i in range(3) if i not in selected_axes][0]
    center_z = (mesh.bounds[iz*2] + mesh.bounds[(iz*2)+1]) / 2
    center_z = 0
    return (xmax + offset, center_y - 0, center_z)


def get_mixing_data(atom_site_info):
    label_mix_map = {}
    loop_vals = [[], [], [], [], [], [], [], []]
    for k, v in atom_site_info.items():
        loop_vals[0].append(k)
        loop_vals[1].append(v['element'])
        loop_vals[2].append(v['symmetry_multiplicity'])
        loop_vals[3].append(v['wyckoff_symbol'])
        loop_vals[4].append(v['x_frac_coord'])
        loop_vals[5].append(v['y_frac_coord'])
        loop_vals[6].append(v['z_frac_coord'])
        loop_vals[7].append(v['site_occupancy'])

    sites = {}
    all_coords = [] 

    for i in range(len(loop_vals[0])):
        all_coords.append([CFloat(loop_vals[j][i]).n for j in range(4, 7)])
    all_coords = np.array(all_coords)

    num_sites = all_coords.shape[0]
    for i in range(num_sites):
        ind_same_coord = (np.linalg.norm(all_coords - all_coords[i, :], axis=1) == 0).astype(int)
        
        if sum(ind_same_coord) == 1:
            if float(loop_vals[-1][i]) == 1.0:
                sites[loop_vals[0][i]] = [loop_vals[j][i] for j in range(1, len(loop_vals))]
            else:
                site_name = loop_vals[0][i]
                label_mix_map[site_name] = site_name
                val = [[loop_vals[1][i], 'Vac'],
                       loop_vals[2][i],
                       loop_vals[3][i],
                       loop_vals[4][i], loop_vals[5][i], loop_vals[6][i],
                       [loop_vals[7][i], 1.0-loop_vals[7][i]]
                       ]
                sites[site_name] = val
        else:
            site_names = [loop_vals[0][j] for j in range(num_sites) if ind_same_coord[j]]
            site_name = '/'.join(site_names)
            # label_mix_map[site_name] = [_site_name for _site_name in site_names]
            for label in site_names:
                label_mix_map[label] = site_name

            if not site_name in sites:
                val = [[loop_vals[1][j] for j in range(num_sites) if ind_same_coord[j]],
                       loop_vals[2][i],
                       loop_vals[3][i],
                       loop_vals[4][i], loop_vals[5][i], loop_vals[6][i],
                       [loop_vals[7][j] for j in range(num_sites) if ind_same_coord[j]]
                       ]
                sites[site_name] = val

    return sites, label_mix_map


def get_colors_for_disorder(sphere, point, mask, colors_fractions, axis_vertical):

    colors_fractions = sorted(colors_fractions, key=lambda x: x[-1], reverse=True)
    
    pts = sphere.points - np.array(point)
    view_dir = (~mask).astype(float)
    
    # Define 12 o'clock as world Z projected onto screen plane
    world_up = np.array([0, 0, 0])
    world_up[axis_vertical] = 1
    
    # 12 o'clock axis: project world_up onto plane perpendicular to view_dir
    up_screen = world_up - np.dot(world_up, view_dir) * view_dir
    up_screen /= np.linalg.norm(up_screen)  # points to 12 o'clock

    # 3 o'clock axis: perpendicular to both view_dir and up_screen
    right_screen = np.cross(up_screen, view_dir)
    right_screen /= np.linalg.norm(right_screen)  # points to 3 o'clock

    # Project points onto screen plane
    u = pts @ right_screen  # horizontal (3 o'clock = positive)
    v = pts @ up_screen     # vertical   (12 o'clock = positive)

    # Angle from 12 o'clock, clockwise
    angles = (np.degrees(np.arctan2(u, v))) % 360

    # Assign wedges
    colors_array = np.zeros(sphere.n_points, dtype=int)
    cumulative = 0.0
    for i, (_, fraction) in enumerate(colors_fractions):
        start_angle = cumulative * 360
        end_angle   = (cumulative + fraction) * 360
        if i == len(colors_fractions) - 1:
            mask = angles >= start_angle
        else:
            mask = (angles >= start_angle) & (angles < end_angle)
        colors_array[mask] = i
        cumulative += fraction
    
    return colors_array


def plot_supercell_pyvista(cif_path, no_label=False, ncols=2, rscale=0.3, fontsize=100, cam_dist=50, cam_tilt=5,
                           width=10000, height=10000, ambient=0.4, diffuse=0.9,
                           theta_resolution=100, phi_resolution=100, manual_cam_pos=None):

    cif = Cif(cif_path)
    lengths = cif.unitcell_lengths
    angles = cif.unitcell_angles
    cif.compute_connections()

    supercell_points_z = cif.supercell_points
    unitcell_points_z = cif.unitcell_points_for_plotting

    unitcell_points_z = [[p[:3], p[3]] for p in unitcell_points_z]
    supercell_points_z = [[p[:3], p[3]] for p in supercell_points_z]

    # convert to cartesian coordinates
    supercell_points = []
    for i in range(len(supercell_points_z)):
        point, s = supercell_points_z[i]
        cart = unit.fractional_to_cartesian(point,
                                            lengths,
                                            angles)
        supercell_points.append((np.array(cart), s))

    unitcell_points = []
    for i in range(len(unitcell_points_z)):
        point, s = unitcell_points_z[i]

        cart = unit.fractional_to_cartesian(point,
                                            lengths,
                                            angles)
        unitcell_points.append((np.array(cart), s))

    loop_vals = cif._loop_values
    site_symbol_map = dict(zip([l for l in loop_vals[0]], [s for s in loop_vals[1]]))

    site_data, label_mix_map = get_mixing_data(cif.atom_site_info)
    
    # Trigger connections computation (lazy property)
    _ = cif.shortest_distance

    selected_CNs = None
    if os.path.isfile("CN.txt"):
        user_prefs = pd.read_csv("CN.txt")
        selected_CNs = dict(zip(user_prefs['Site'].tolist(), user_prefs['CN'].tolist()))
        print(selected_CNs)
        print(f"CN values from CN.txt will be used for sites lited in CN.txt")

    colors = {el: get_element_color(el) for el in set(list(site_symbol_map.values()))}
    colors['Vac'] = 'white'
    if os.path.isfile("colors.txt"):
        user_prefs = pd.read_csv("colors.txt")
        slected_colors = dict(zip(user_prefs['Site'].tolist(), user_prefs['Color'].tolist()))
        for k in colors:
            if k in slected_colors:
                colors[k] = slected_colors[k]
    
    plotter = pv.Plotter(window_size=(width, height))


    # setup mask for 2D
    mask = np.array([False, False, False])
    sorted_indices = np.argsort(lengths)[::-1]
    largest_indices = sorted_indices[:2]
    mask[largest_indices] = True

    selected_axes = [i for i in range(3) if mask[i]]
    axis_vertical, axis_horizontal = selected_axes
    if lengths[selected_axes[0]] >= lengths[selected_axes[1]]:
        axis_horizontal, axis_vertical = selected_axes

    camera_orientation = [0, 0, 0]
    camera_pos = [0, 0, 0]
    camera_pos[np.argwhere(mask == False).squeeze()] = cam_dist 
    camera_orientation[axis_vertical] = 1
    plotter.camera.up = camera_orientation
    plotter.camera.position = camera_pos

    plotter.enable_parallel_projection()
    plotter.reset_camera()

    # Plot each atom as a sphere
    for point, label in unitcell_points:
        
        if label in label_mix_map:
            mlabel = label_mix_map[label]
            msite = site_data[mlabel]
            colors_fractions = [(elem, occ) for elem, occ in zip(msite[0], msite[-1])]
            colors_fractions = sorted(colors_fractions, key=lambda x: element_data[x[0]][0])
            colors_fractions = [(colors.get(v[0], 'black'), v[1]) for v in colors_fractions]

            element = msite[0][0]
            if msite[-1][0] < msite[-1][1]:
                element = msite[0][1]

            sphere = pv.Sphere(center=point, radius=element_data[element][1]*rscale, theta_resolution=theta_resolution, 
                               phi_resolution=phi_resolution)
            colors_array = get_colors_for_disorder(sphere, point, mask, colors_fractions, axis_vertical)
            sphere["region"] = colors_array.astype(float)
            cmap = ListedColormap([c for c, _ in colors_fractions])
            plotter.add_mesh(sphere, scalars="region", cmap=cmap,
                        clim=[0, len(colors_fractions) - 1],
                        show_scalar_bar=False, smooth_shading=True, ambient=ambient, diffuse=diffuse)
        else:
            element = site_symbol_map.get(label, '?')
            color = colors.get(element, 'black')
            sphere = pv.Sphere(center=point, radius=element_data[element][1]*rscale)
            plotter.add_mesh(sphere, color=color, show_scalar_bar=True, ambient=ambient, diffuse=diffuse)

    # add box
    unitcell_hull = [
        [0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
        [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]
    ]

    unitcell_hull_cartesian = [unit.fractional_to_cartesian(point, lengths, angles)
                               for point in unitcell_hull]

    edges = [(0, 2), (0, 4), (2, 6), (4, 6),
             (0, 1), (1, 5), (4, 5), 
             (1, 3), (2, 3), 
             (3, 7),
             (4, 6), 
             (5, 7), 
             (6, 7)
             ]
    
    for i, j in edges:
        pi = unitcell_hull_cartesian[i]
        pj = unitcell_hull_cartesian[j]
        cylinder = pv.Cylinder(center=(pi+pj)/2, direction=pj-pi, height=np.linalg.norm(pj-pi), radius=0.01, resolution=16)
        plotter.add_mesh(cylinder, color='black')

    axes = {'x': axis_horizontal, 'y': axis_vertical}

    a, b, c = lengths
    alpha, beta, gamma = angles

    # Build the fractional -> Cartesian matrix (standard crystallographic convention)
    cos_a, cos_b, cos_g = np.cos(alpha), np.cos(beta), np.cos(gamma)
    sin_g = np.sin(gamma)

    vol_factor = np.sqrt(
        1 - cos_a**2 - cos_b**2 - cos_g**2
        + 2 * cos_a * cos_b * cos_g
    )

    M = np.array([
        [a,          b * cos_g,   c * cos_b],
        [0,          b * sin_g,   c * (cos_a - cos_b * cos_g) / sin_g],
        [0,          0,           c * vol_factor / sin_g]
    ])

    conns = cif.connections
 
    for i, site in enumerate(site_data.keys(), 1):
        site_label = site
        if site in label_mix_map.values():
            site = [k for k in label_mix_map.keys() if label_mix_map[k] == site][0]

        points_wd = conns[site][:21]
        
        if selected_CNs and site_label in selected_CNs:
            CN = selected_CNs[site_label]
        else:
            CN = CN_of_site(points_wd)

        row = int(i / ncols)
        col = i % ncols

        translation_vector = [0, 0, 0]
        translation_vector[axes['x']] = col
        translation_vector[axes['y']] = -row

        t_cart = M @ np.array(translation_vector)

        _neighbors = [[v[-1], v[0]] for v in points_wd[:CN]]

        # translate
        neighbors = [[np.array(points_wd[0][2])+t_cart, site]]
        for val in _neighbors:
            val[0] = np.array(val[0]) + t_cart
            neighbors.append(val)

        
        for point, label in neighbors:
            if label in label_mix_map:

                mlabel = label_mix_map[label]
                msite = site_data[mlabel]

                colors_fractions = [(colors.get(elem, 'black'), occ) for elem, occ in zip(msite[0], msite[-1])]

                elements = msite[0]
                sphere = pv.Sphere(center=point, radius=element_data[elements[0]][1]*rscale, theta_resolution=theta_resolution, 
                               phi_resolution=phi_resolution)

                colors_array = get_colors_for_disorder(sphere, point, mask, colors_fractions, axis_vertical)
                sphere["region"] = colors_array.astype(float)
                cmap = ListedColormap([c for c, _ in colors_fractions])
                plotter.add_mesh(sphere, scalars="region", cmap=cmap,
                            clim=[0, len(colors_fractions) - 1],
                            show_scalar_bar=False, smooth_shading=True, ambient=ambient, diffuse=diffuse)
                
                label = mlabel
            else:
                element = site_symbol_map.get(label, '?')
                color = colors.get(element, 'black')

                sphere = pv.Sphere(center=point, radius=element_data[element][1]*rscale)
                plotter.add_mesh(sphere, color=color, show_scalar_bar=True, ambient=ambient, diffuse=diffuse)


        points = np.array([neighbor[0] for neighbor in neighbors])
        hull = ConvexHull(points)
        
        # Build mesh from hull triangles (2D faces of the 3D hull)
        faces = []
        for simplex in hull.simplices:
            face = np.array([3, simplex[0], simplex[1], simplex[2]], dtype=np.int64)
            faces.append(face)
        
        faces = np.concatenate(faces)
        polyhedron = pv.PolyData(points, faces=faces)

        if site in label_mix_map:
            msite = site_data[label_mix_map[site]]
            element = msite[0][0]
            if msite[-1][0] < msite[-1][1]:
                element = msite[0][1]
        else:
            element = site_symbol_map[site]

        plotter.add_mesh(polyhedron, color=colors[element], opacity=0.3, show_edges=False, ambient=ambient, diffuse=diffuse)
        edges = polyhedron.extract_feature_edges()
        plotter.add_mesh(edges, color="black", line_width=2, opacity=0.7)

        if not no_label:
            label_coord = neighbors[0][0].copy() 
            min_x = np.min(points[:, axis_horizontal])
            max_x = np.max(points[:, axis_horizontal])
            min_y = np.min(points[:, axis_vertical])
            max_y = np.max(points[:, axis_vertical])
            
            # Position: Top-Right corner of the bounding box, with small padding
            label_coord[axis_horizontal] = max_x - (abs(max_x-min_x)*0.9) # 0.2 
            label_coord[axis_vertical] = max_y + (abs(max_y-min_y)*0.05) # 0.2

            plotter.add_point_labels(label_coord, [site_label],
                                    font_size=fontsize,
                                    text_color="black",
                                    bold=False,
                                    show_points=False,
                                    always_visible=True,
                                    shape=None,)
        # break
        
    # element legends
    unitcell_hull_cartesian = np.array(unitcell_hull_cartesian)
    ax_out = [i for i in range(3) if i not in [axis_horizontal, axis_vertical]][0]

    left_bottom_mid = [0, 0, 0]
    left_bottom_mid[ax_out] = 0.5
    left_bottom_mid = unit.fractional_to_cartesian(left_bottom_mid, lengths, angles)

    left_top_mid = [0, 0, 0]
    left_top_mid[ax_out] = 0.5
    left_top_mid[axis_vertical] = 1.0
    left_top_mid = unit.fractional_to_cartesian(left_top_mid, lengths, angles)

    anchor_point = left_top_mid
    if left_bottom_mid[axis_horizontal] < left_top_mid[axis_horizontal]:
        anchor_point = left_bottom_mid
        left_top_mid[axis_horizontal] = left_bottom_mid[axis_horizontal]

    colors = dict(sorted(colors.items(), key=lambda item: element_data[item[0]][2], reverse=True))
    for i, (el, color) in enumerate(colors.items()):
        if el == 'Vac': continue
        center = anchor_point.copy()
        center[axis_horizontal] -= lengths[axis_horizontal] * 0.3
        center[axis_vertical] -= (lengths[axis_vertical] * 0.5 + i*1.75) 
        sphere = pv.Sphere(center=center, radius=element_data[element][1]*rscale)
        plotter.add_mesh(sphere, color=color, show_scalar_bar=True, ambient=ambient, diffuse=diffuse)
        
        if not no_label:
            center[axis_horizontal] *= 1.5
            center[axis_vertical] *= 0.9

            plotter.add_point_labels(center, [el],
                                    font_size=fontsize,
                                    text_color="black",
                                    bold=False,
                                    show_points=False,
                                    always_visible=True,
                                    shape=None,)
    
    axis_labels = ['a', 'b', 'c']
    left_top_mid[axis_horizontal] -= lengths[axis_horizontal] * 0.4
    left_top_mid[axis_vertical] -= lengths[axis_vertical] * 0.3
    
    for i, (color, ax) in enumerate(zip(['red', 'green', 'blue'], [axis_horizontal, axis_vertical, ax_out])):
        direction = np.array([0, 0, 0])
        direction[i] = 1
        arrow = pv.Arrow(start=left_top_mid, direction=direction, scale=2, 
                         tip_length=0.3, tip_radius=0.1, shaft_radius=0.05,)
        plotter.add_mesh(arrow, color=color)

        if not no_label:
            plotter.add_point_labels(left_top_mid+(direction*2.2), [axis_labels[ax]],
                                        font_size=fontsize,
                                        text_color="black",
                                        bold=False,
                                        show_points=False,
                                        always_visible=True,
                                        shape=None,)


    plotter.reset_camera(bounds=None)   # auto-fit
    plotter.camera.zoom(0.9)

    # Save the plot
    filename = f"unitcell_and_sites_plot"
    if no_label:
        filename += "_nolabel"
    plotter.export_html("scene.html")
    plotter.screenshot(f'{filename}.png')
    print(f"Screenshot saved to {filename}.png")

    camera_pos = np.array(camera_pos)
    
    if manual_cam_pos is None:
        camera_pos[camera_pos==0] = cam_tilt
    else:
        camera_pos = manual_cam_pos

    plotter.reset_camera()
    plotter.camera.position = camera_pos
    plotter.reset_camera()
    print(f"Tilt position used: {camera_pos}")

    plotter.screenshot(f'{filename}_tilt.png')
    print(f"Screenshot saved to {filename}_tilt.png")
    plotter.close()


def cli_plot_supercell_pyvista():
    """CLI entry point for plotting supercell structures."""
    parser = argparse.ArgumentParser(description="Plot crystal structure using PyVista.")
    
    parser.add_argument("cif_path", help="Path to the CIF file")
    parser.add_argument("--nolabels", action='store_true', help="Generate images without labels")
    parser.add_argument("--ncols", type=int, default=2, help="Number of columns in the plot grid")
    parser.add_argument("--rscale", type=float, default=0.3, help="Scale factor for atom radii")
    parser.add_argument("--fontsize", type=int, default=100, help="Font size for labels")
    parser.add_argument("--cam-dist", type=float, default=50, help="Camera distance from origin")
    parser.add_argument("--cam-tilt", type=float, default=5, help="Camera tilt angle")
    parser.add_argument("--width", type=int, default=8000, help="Image width in pixels")
    parser.add_argument("--height", type=int, default=8000, help="Image height in pixels")
    parser.add_argument("--cam-pos", type=lambda s: ast.literal_eval(s), default=None, help="Camera position e.g. '[100, 5, 5]' for a structure viewed through [001] with small tilt in a and b axis")
    
    args = parser.parse_args()
    
    plot_supercell_pyvista(
        cif_path=args.cif_path,
        no_label=args.nolabels,
        ncols=args.ncols,
        rscale=args.rscale,
        fontsize=args.fontsize,
        cam_dist=args.cam_dist,
        cam_tilt=args.cam_tilt,
        width=args.width,
        height=args.height,
        manual_cam_pos=args.cam_pos
    )


if __name__ == "__main__":
    # cif_path = "/home/bala/Documents/44_shelxl/test_files/2_leaves_done/test.cif"
    cif_path = "/home/bala/Documents/44_shelxl/test_files/4_DyIrSn_done/test.cif"
    # cif_path = "1232634.cif"
    plot_supercell_pyvista(cif_path, cam_tilt=8, fontsize=60, width=3000, height=3000)
