import json
import math
import random
from copy import deepcopy
import numpy as np

AVOGADRO = 6.02214076e23

def build_system(
        bond_list, 
        # dihedrals,
        xyz_file, 
        n_molecule, 
        # density,
        outfile,
        ):
    """
        units in:
            angstrom
        units json out:
            nm
        units xyz out:
            angstrom
        units internally
            angstrom 
    """
    typelist = {"O": "O-TIP3P", "H": "H-TIP3P", "Ar": "Ar"}

    elements, coords_angst = load_xyz(xyz_file)
    n_atoms = len(elements)
    bonds = []
    atoms = []
    box_size_angst = np.array([20, 20, 20])
    molecule_dimensions_angst, centered_coords_angst = center_molecule(coords=coords_angst)
    print(centered_coords_angst)

    total_coords = np.zeros((n_molecule * n_atoms,3))
    total_elements = []

    for i in range(n_molecule):
        total_elements += elements 

        rotated_molecule = np.einsum('ij, nj-> ni', random_rotation_matrix(), centered_coords_angst)

        shiftx = random.uniform(-(0.5*box_size_angst[0] - 2 - molecule_dimensions_angst[0]/2), (0.5*box_size_angst[0] - 2 - molecule_dimensions_angst[0]/2))
        shifty = random.uniform(-(0.5*box_size_angst[1] - 2 - molecule_dimensions_angst[1]/2), (0.5*box_size_angst[1] - 2 - molecule_dimensions_angst[1]/2))
        shiftz = random.uniform(-(0.5*box_size_angst[2] - 2 - molecule_dimensions_angst[2]/2), (0.5*box_size_angst[2] - 2 - molecule_dimensions_angst[2]/2))
        random_shift = [shiftx, shifty, shiftz]
        total_coords[i*n_atoms:(i+1)*n_atoms] = shift_molecule(rotated_molecule, random_shift)

        new_bonds = np.array(bond_list) + i*n_atoms
        new_bonds = [[index + i * n_atoms for index in row] for row in bond_list]
        [bonds.append(b) for b in new_bonds]

    print(bonds)

    write_xyz(elements=total_elements, coords=total_coords)

    data = {"box": list(angstrom_to_nm(box_size_angst))}
    for i,elem in enumerate(total_elements):
        print(i)
        atom = {
                "id": i, 
                "type": typelist[elem], 
                "pos": list(angstrom_to_nm(total_coords[i]))
                }
        atoms.append(atom)
    data["atoms"] = atoms
    data["bonds"] = bonds

    with open(f"{outfile}.json", "w", encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)
    

def load_xyz(filepath):
    with open(filepath, "r") as f:
        line1 = f.readline()
        N = int(line1)
        coords = np.zeros(shape=(N,3))
        next(f)
        atoms = []

        for i,line in enumerate(f):
            element, x, y, z = line.split()
            element = str(element)
            x = float(x)
            y = float(y)
            z = float(z)
            atoms.append(element)
            coords[i] = [x,y,z]
    assert len(atoms) == coords.shape[0]
    return atoms, coords

def shift_molecule(coords_in, shift_vec):
    return coords_in + shift_vec

def center_molecule(coords):
    min_xyz = np.zeros(3)
    len_xyz = np.zeros(3)
    for i in range(3):
        min_xyz[i] = np.min(coords[:,i])
        len_xyz[i] = np.max(coords[:,i]) - np.min(coords[:,i])
    center_of_volume = min_xyz - len_xyz/2 
    return len_xyz, shift_molecule(coords_in=coords, shift_vec=-center_of_volume)
    
def write_xyz(elements, coords):
    with open("test.xyz", "w") as f:
        f.write(f"{len(elements)}\n")
        f.write("comment_line\n")
        for i,elem in enumerate(elements):
            f.write(f"{elem}\t{coords[i,0]}\t{coords[i,1]}\t{coords[i,2]}\n")

def angstrom_to_nm(angstrom):
    return 0.1 * angstrom
        
def random_rotation_matrix():
    """
    Generate a random 3D rotation matrix.
    """

    u1 = random.random()
    u2 = random.random()
    u3 = random.random()

    q1 = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    q2 = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2 * math.pi * u3)

    R = np.array([
        [1 - 2*(q3*q3 + q4*q4), 2*(q2*q3 - q1*q4),     2*(q2*q4 + q1*q3)],
        [2*(q2*q3 + q1*q4),     1 - 2*(q2*q2 + q4*q4), 2*(q3*q4 - q1*q2)],
        [2*(q2*q4 - q1*q3),     2*(q3*q4 + q1*q2),     1 - 2*(q2*q2 + q3*q3)]
    ])

    return R

def read_types():
    pass



if __name__ == "__main__":
    bond_list = [[0,1], [0,2]]
    # bond_list = []
    build_system(outfile="test", bond_list=bond_list, xyz_file="h2o.xyz", n_molecule=20) 