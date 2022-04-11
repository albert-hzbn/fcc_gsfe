"""
This program reads FCC structure in the POSCAR format in cartisian coordinates,
introduces stacking fault in them and output POSCAR file.
"""

import math

len_stacking = int(input("Enter number of layers: "))

input_range = input(
    "Enter range of top layer you want ot shift for isf (e.g. entering 6-9 shifts layer from 6 to 9): ")
s_i, eI = input_range.split('-')
s_i, eI = int(s_i)-1, int(eI)

max_shift = int(
    input("Enter maximum n(whole number) to shift by magnitude n*b_twin: "))


# Reading ideal disordered POSCAR file
data = []
with open('POSCAR', 'r') as reader:
    data = reader.read().split("\n")

# Extracting elements, number of atoms and position
elem_arr = data[5].split()
num_atoms_arr = list(map(lambda x: int(x), data[6].split()))
atoms = sum(num_atoms_arr)
atoms_str = data[8:]
atoms_arr = []
elem_arr_len = len(elem_arr)
i = 0
for e_i in range(0, elem_arr_len):
    for a_i in range(0, num_atoms_arr[e_i]):
        atom_pos = list(map(lambda x: float(x), atoms_str[i].split()))
        atoms_arr.append([elem_arr[e_i]]+atom_pos)
        i += 1

# Sorting data
atoms_arr = sorted(atoms_arr, key=lambda x: x[3])

# Lattive vectors
a1, a2, a3 = list(map(lambda x: float(x), data[2].split()))
b1, b2, b3 = list(map(lambda x: float(x), data[3].split()))
c1, c2, c3 = list(map(lambda x: float(x), data[4].split()))


# Converts data to POSCAR
def convert_to_poscar(data, atoms_arr, current_shift):

    shift_x = current_shift*math.cos(30*math.pi/180)
    shift_y = current_shift*math.sin(30*math.pi/180)

    # For shifting lattice vectors
    c1, c2 = shift_x, shift_y

    symbol_str = ' '.join([str(elem) for elem in elem_arr])
    num_atoms_str = ' '.join([str(n) for n in num_atoms_arr])

    # Adding lattice vectors and some other info. to the array
    new_data = [
        data[0],
        data[1],
        f'  {a1:.10f}  {a2:.10f}  {a3:.10f}',
        f'  {b1:.10f}  {b2:.10f}  {b3:.10f}',
        f'  {c1:.10f}  {c2:.10f}  {c3:.10f}',
        f'  {symbol_str}',
        f'  {num_atoms_str}',
        "Selective Dynamics",
        "Cartesian"
    ]

    new_atoms_arr = []
    for e_i in range(0, elem_arr_len):
        for a_i in range(0, atoms):
            if(elem_arr[e_i] == atoms_arr[a_i][0]):
                new_atoms_arr.append(atoms_arr[a_i])

    # For atomic positions
    for a_i in range(0, atoms):
        new_data.append(
            f'  {new_atoms_arr[a_i][1]:.10f}  {new_atoms_arr[a_i][2]:.10f}  {new_atoms_arr[a_i][3]:.10f}  F F T')
    return new_data


# Saves data to file
def save_to_file(f_name, data):
    with open(f_name, 'w') as f:
        f.write('\n'.join(data))


# Magnitude of burger vector
b_twin = 1/math.sqrt(6)

# minimum  shift
shift = 0.1*b_twin

# Minimum shift in x and y direction
shift_x = shift*math.cos(30*math.pi/180)
shift_y = shift*math.sin(30*math.pi/180)

num_atoms_layer = int(atoms/len_stacking)

# Loop for shifting atoms
for ind in range(0, max_shift):
    for f in range(0, 10):
        for layer in range(0, len_stacking):

            # Condition for range of layers to be shifted
            if((layer+1) >= (s_i+ind) and (layer+1) <= eI):
                li = num_atoms_layer*layer
                ri = num_atoms_layer*(layer+1)
                for i in range(li, ri):
                    atoms_arr[i][1] += shift_x
                    atoms_arr[i][2] += shift_y

        current_shift = ((f+1)+ind*10)/10*b_twin
        new_data = convert_to_poscar(data, atoms_arr, current_shift)
        save_to_file("POSCAR_"+str((f+1)+ind*10), new_data)
