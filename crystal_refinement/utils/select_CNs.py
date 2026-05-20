from cifkit import Cif
from crystal_refinement.utils.site_plots import CN_of_site, get_mixing_data, get_element_color
import argparse

web_colors = [
    'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 
    'black', 'blanchedalmond', 'blue', 'blueviolet', 'brown', 'burlywood', 
    'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 
    'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 
    'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 
    'darkorange', 'darkorchid', 'darkred', 'darksalmon', 'darkseagreen', 
    'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 
    'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 
    'firebrick', 'floralwhite', 'forestgreen', 'fuchsia', 'gainsboro', 
    'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 
    'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 
    'lavenderblush', 'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 
    'lightcyan', 'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey', 
    'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 
    'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 
    'linen', 'magenta', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid', 
    'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen', 
    'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 
    'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 
    'orangered', 'orchid', 'palegoldenrod', 'palegreen', 'paleturquoise', 
    'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 
    'powderblue', 'purple', 'rebeccapurple', 'red', 'rosybrown', 'royalblue', 
    'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 
    'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'snow', 
    'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 
    'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen'
]


def select_CNs(cif_path=None, with_colors=False):
    # Handle CLI arguments if called directly
    if cif_path is None:
        parser = argparse.ArgumentParser(description="Select CNs and optionally colors for crystal sites")
        parser.add_argument("cif_path", help="Path to the CIF file")
        parser.add_argument("-c", "--colors", action="store_true", help="Enable color selection for sites")
        args = parser.parse_args()
        cif_path = args.cif_path
        with_colors = args.colors

    cif = Cif(cif_path)
    cif.compute_connections()

    loop_vals = cif._loop_values
    site_symbol_map = dict(zip([l for l in loop_vals[0]], [s for s in loop_vals[1]]))

    conns = cif.connections
    
    site_data, label_mix_map = get_mixing_data(cif.atom_site_info)
    # print(site_data) # Optional: debug info
    # print(label_mix_map) # Optional: debug info

    # Store results to display all at once
    site_cns = {}
    for site, _ in site_data.items():
        site_label = site
        if site in label_mix_map.values():
            site_label = [k for k in label_mix_map if label_mix_map[k]==site][0]

        neighbors = conns[site_label]
        site_cns[site] = CN_of_site(neighbors[:21], full=True)

    # ---------------------------------------------------------
    # PHASE 1: Coordination Number Selection
    # ---------------------------------------------------------
    print("=== Phase 1: Select Coordination Numbers ===")
    print("-" * 40)
    print(f"{'Site':<10}: CN Values")
    print("-" * 40)
    for site, cn_vals in sorted(site_cns.items()):
        print(f"{site:<10}: {cn_vals}")
    print("-" * 40)
    
    print("\nEnter a CN for each site.")
    print("If you press Enter without typing, the first value will be used.\n")

    # Prompt user for CN selection only
    selections = {}
    
    for site in sorted(site_cns.keys()):
        cn_vals = site_cns[site]
        default_cn = cn_vals[0]
        
        selected_cn = None
        max_tries = 3

        if site in label_mix_map.values():
            site_label = [k for k in label_mix_map if label_mix_map[k]==site][0]
        else:
            site_label = site
        
        for attempt in range(1, max_tries + 1):
            prompt = f"Enter CN for {site} (default: {default_cn}): "
            if attempt > 1:
                prompt += f" (Attempt {attempt}/{max_tries}, valid options: {cn_vals}) "
            
            user_input = input(prompt).strip()
            
            if user_input == "":
                selected_cn = default_cn
                break
            
            try:
                candidate_cn = int(user_input)
                if candidate_cn in cn_vals:
                    selected_cn = candidate_cn
                    break
                else:
                    print(f"Error: {candidate_cn} is not in the computed CNs for {site}.")
                    if attempt == max_tries:
                        print(f"Max attempts reached. Using default value: {default_cn}")
                        selected_cn = default_cn
            except ValueError:
                print(f"Error: '{user_input}' is not a valid integer.")
                if attempt == max_tries:
                    print(f"Max attempts reached. Using default value: {default_cn}")
                    selected_cn = default_cn
        
        selections[site] = selected_cn

    # Write CNs to CN.txt
    cn_output_file = "CN.txt"
    with open(cn_output_file, 'w') as f:
        f.write(f"Site,CN\n")
        for site in sorted(selections.keys()):
            cn = selections[site]
            f.write(f"{site},{cn}\n")
    
    print(f"\nCoordination numbers written to {cn_output_file}")
    with open(cn_output_file, 'r') as f:
        print("Contents:")
        print(f.read())

    # ---------------------------------------------------------
    # PHASE 2: Color Selection (if enabled)
    # ---------------------------------------------------------
    if with_colors:
        print("\n" + "=" * 40)
        print("=== Phase 2: Select Colors ===")
        print("-" * 40)
        
        print("\nWeb colors available for selection:")
        print(", ".join(web_colors[:10]))
        print("... (and more)")
        print("See https://www.w3.org/TR/css-color-4/#named-colors")
        print(f"Total: {len(web_colors)} colors\n")

        color_selections = {}

        # Get unique elements from site_symbol_map to avoid duplicate prompts for same element
        # We map Element -> List of Sites having that element
        element_sites = {}
        for site, element in site_symbol_map.items():
            if element not in element_sites:
                element_sites[element] = []
            element_sites[element].append(site)

        print(f"Selecting colors for {len(element_sites)} unique element(s).\n")

        for element, sites in sorted(element_sites.items()):
            default_color = get_element_color(element)
            selected_color = None
            max_color_tries = 3
            
            print(f"Element '{element}' is found in sites: {', '.join(sites)}")

            for attempt in range(1, max_color_tries + 1):   
                prompt = f"Enter color for element {element} (default: {default_color}): "
                if attempt > 1:
                    prompt += f" (Attempt {attempt}/{max_color_tries}, valid options: any of {len(web_colors)} web colors) "
                
                user_input = input(prompt).strip()
                
                if user_input == "":
                    selected_color = default_color
                    break
                
                # Normalize input for comparison
                normalized_input = user_input.lower().replace(" ", "")
                
                if normalized_input in [c.replace(" ", "").lower() for c in web_colors]:
                    # Find the original case-sensitive version
                    original_color = next(c for c in web_colors if c.replace(" ", "").lower() == normalized_input)
                    selected_color = original_color
                    break
                else:
                    print(f"Error: '{user_input}' is not a valid web color.")
                    if attempt == max_color_tries:
                        print(f"Max attempts reached. Using default color: {default_color}")
                        selected_color = default_color
            
            # Assign this color to ALL sites that have this element
            # for site in sites:
            color_selections[element] = selected_color

        # Write Colors to colors.txt
        colors_output_file = "colors.txt"
        with open(colors_output_file, 'w') as f:
            f.write(f"Site,Color\n")
            for site in sorted(color_selections.keys()):
                color = color_selections[site]
                f.write(f"{site},{color}\n")
        
        print(f"\nColors written to {colors_output_file}")
        with open(colors_output_file, 'r') as f:
            print("Contents:")
            print(f.read())


if __name__ == "__main__":
    select_CNs()