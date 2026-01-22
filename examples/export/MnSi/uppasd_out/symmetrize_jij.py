import numpy as np


def symmetrize_jij(input_file, output_file, tol=1e-6):
    """
    Reads Jij data, preserves all entries, and for entries where the
    position vector r_ij is (within tol) the negative of another entry's
    r_ji (and indices are swapped), averages the Jij values across those
    matched entries. Each original entry keeps its r_x/r_y/r_z and r_ij.
    """

    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    if not lines:
        print(f"Input file {input_file} is empty")
        return

    # Parse all lines into a list; preserve malformed or empty lines
    jij_list = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            jij_list.append(None)
            continue

        tokens = stripped.split()
        if len(tokens) < 7:
            jij_list.append({'raw': line})
            continue

        try:
            i = int(tokens[0])
            j = int(tokens[1])
            r_x = float(tokens[2])
            r_y = float(tokens[3])
            r_z = float(tokens[4])
            jij = float(tokens[5])
            r_ij = float(tokens[6])
        except ValueError:
            jij_list.append({'raw': line})
            continue

        jij_list.append({
            'i': i, 'j': j, 'r_x': r_x, 'r_y': r_y, 'r_z': r_z,
            'jij': jij, 'r_ij': r_ij
        })

    n = len(jij_list)
    processed = [False] * n
    jij_new = [None] * n

    for idx, entry in enumerate(jij_list):
        if processed[idx] or not entry or 'raw' in entry:
            continue

        i = entry['i']
        j = entry['j']
        rx = entry['r_x']
        ry = entry['r_y']
        rz = entry['r_z']

        # gather matches: later entries with swapped indices and negated r vector
        matches = [idx]
        for k in range(idx + 1, n):
            other = jij_list[k]
            if not other or 'raw' in other or processed[k]:
                continue
            if other.get('i') == j and other.get('j') == i:
                if np.allclose([rx, ry, rz], [-other['r_x'], -other['r_y'], -other['r_z']], atol=tol):
                    matches.append(k)

        vals = [jij_list[m]['jij'] for m in matches]
        avg = float(np.mean(vals)) if vals else 0.0
        for m in matches:
            jij_new[m] = avg
            processed[m] = True

    # remaining data entries that had no matches keep their original Jij
    for idx, entry in enumerate(jij_list):
        if entry and 'raw' not in entry and not processed[idx]:
            jij_new[idx] = entry['jij']
            processed[idx] = True

    # Write output preserving original order
    with open(output_file, 'w', encoding='utf-8') as out:
        for idx, entry in enumerate(jij_list):
            if entry is None:
                out.write('\n')
            elif 'raw' in entry:
                out.write(entry['raw'])
            else:
                i = entry['i']
                j = entry['j']
                rx = entry['r_x']
                ry = entry['r_y']
                rz = entry['r_z']
                r_ij = entry['r_ij']
                val = jij_new[idx]
                out.write(f"{i} {j} {rx:.6f} {ry:.6f} {rz:.6f} {val:.6f} {r_ij:.6f}\n")

    print(f"Symmetrized Jij interactions saved to {output_file}")


if __name__ == '__main__':
    symmetrize_jij('jfile', 'jfile.symmed')
