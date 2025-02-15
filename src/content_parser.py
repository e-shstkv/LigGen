import math
import pKa
import config as cfg

"""
moduleContent parser - get info from .cif file format.
"""


def contents_pdb(cif_file_path=None):
    """Formation lines from files"""
    if cif_file_path:
        with open(cif_file_path, 'r') as r:
            lines = []
            content = r.readlines()
            for line in content:
                line = line.split()
                lines.append(line)

            return lines
    else:
        print(f"pdb_content: {rcsb}")


def get_pdb_id():
    cif_file_path = cfg.configuration_read('file input')
    content = contents_pdb(cif_file_path)
    answer = 0
    for line in content:
        if line[0] == '_entry.id':
            answer = line[-1]

        elif line[0] == 'PDB':
            answer = line[-1]

        elif line[0] == '_pdbx_database_status.entry_id':
            answer = line[-1]

        elif line[0] == '_cell.entry_id':
            answer = line[-1]

        elif line[0] == '_symmetry.entry_id':
            answer = line[-1]

        elif line[0] == '_struct.entry_id':
            answer = line[-1]

    return answer


def protein_atoms():
    """Preparation "ATOM" lines"""
    cif_file_path = cfg.configuration_read('file input')
    content = contents_pdb(cif_file_path)
    if content:
        atoms = []
        for line in content:
            if line[0] == 'ATOM':
                atom = [list(map(int, line[1:2])), line[2], line[3], line[5], line[8],
                        list(map(float, line[10:13])), ]
                atoms.append(atom)

        return atoms
    else:
        print(f"Protein atoms content: {content}")


def all_hetatm_atoms():
    """Preparation LIST of HETATMS"""
    cif_file_path = cfg.configuration_read('file input')
    content = contents_pdb(cif_file_path)
    hetatms = []

    if content:
        for line in content:
            if line[0] == 'HETATM':
                hetatm = [list(map(int, line[1:2])), line[2], line[3], line[5], line[8],
                          list(map(float, line[10:13])), line[16], ]
                #   If line[8] = true => modified aminoacid by aminoacid; false => true "HETATM"
                hetatms.append(hetatm)
        return hetatms
    else:
        print(f"All hetatm atoms content:{content}")


def h2o():
    """Return only h2o strings"""
    content = all_hetatm_atoms()
    h2os = []
    if content:
        waters = ['HOH', 'DOH', 'HOD', 'DOD']
        for line in content:
            if line[3] in waters:
                h2os.append(line)

        return h2os
    else:
        print(f"Oxygen content:{content}")


def find_empirical_water_env():
    """find practice water environment of aminoacids nearest [angstrom] angstrom to water"""
    angstrom = 5.0
    p_atoms = protein_atoms()
    h2os = h2o()
    pdb_id = get_pdb_id()
    if p_atoms:
        if h2os:
            water_environment_fragment = []
            for oxygen in h2os:
                x0, y0, z0 = oxygen[5][0], oxygen[5][1], oxygen[5][2]
                for p_atom0 in p_atoms:
                    x1, y1, z1 = p_atom0[5][0], p_atom0[5][1], p_atom0[5][2]
                    water_atom_length = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2 + (z0 - z1) ** 2)
                    if water_atom_length < angstrom:
                        if p_atom0 not in water_environment_fragment:
                            water_environment_fragment.append(p_atom0)

            for env in water_environment_fragment:
                for p_atom1 in p_atoms:
                    if env[3] == p_atom1[3]:
                        if env[4] == p_atom1[4]:
                            if p_atom1 not in water_environment_fragment:
                                water_environment_fragment.append(p_atom1)

            # # dump content //:

            # file_name = f"15_02_25_practice_water_env_{pdb_id}.txt"
            # # fip = f"c:\\users\\алексей\\pycharmprojects\\pdbstruct\\3blr.cif"
            # fip = "C:\\Users\\Алексей\\Desktop\\PURE .CIF FILES\\3blr.cif"
            #
            # file_input = contents_pdb(cif_file_path=fip)
            # with open(file_name, 'w') as f:
            #     for water in water_environment_fragment:
            #         # print(water[0], type(water[0]))
            #         for atom in file_input:
            #             if atom[0] == 'ATOM':
            #                 # print()
            #                 # print(water[0][0], atom[1], type(water[0][0]), type(atom[1]))
            #                 if water[0][0] == int(atom[1]):
            #                     # print(water[0][0], atom[1], type(water[0][0]), type(atom[1]))
            #                     atom = f"ATOM   {atom[1]} {atom[2]}  {atom[3]}   {atom[4]} {atom[5]} {atom[6]} {atom[7]} {atom[8]}  ? {atom[10]}  {atom[11]}  {atom[12]} {atom[13]} {atom[14]}  ? {atom[-5]} {atom[-4]} {atom[-3]} {atom[-2]}   {atom[-1]} "
            #                     # atom = f"{atom}"
            #                     f.write(atom + '\n')

            # # dump content \\
            return water_environment_fragment

        else:
            print(f"h2os: {h2os}")
    else:
        print(f"p_atoms:{p_atoms}")


def find_ion_pairs_3_2():

    LG_mode = cfg.configuration_read('mode')
    p_atoms = 0
    if LG_mode == 'prp':
        p_atoms = protein_atoms()
    elif LG_mode == 'prw':
        p_atoms = find_empirical_water_env()
    flag_d = cfg.flag_distant()
    c_a = pKa.conjugate_acids()
    p_atoms_copy = p_atoms[:]
    ion_pairs = []
    """First loop"""
    for i, p_atom0 in enumerate(p_atoms):
        alone_fg = []
        """Second loop"""
        for k0, v0 in c_a.items():
            if p_atom0[3] in k0:
                if p_atom0[2] in v0:
                    for p_atom1 in p_atoms_copy:
                        for k1, v1 in c_a.items():
                            if p_atom1[3] in k1:
                                if p_atom1[2] in v1:
                                    if p_atom1 != p_atom0:
                                        x0, y0, z0 = p_atom0[5][0], p_atom0[5][1], p_atom0[5][2]
                                        x1, y1, z1 = p_atom1[5][0], p_atom1[5][1], p_atom1[5][2]
                                        length0 = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2 + (z0 - z1) ** 2)
                                        if 2.4 < length0 <= 7:
                                            blijniy_4_2_A = 0
                                            dalniy_8_0_A = 0
                                            alone = 0

                                            if p_atom0[3] != 'MET' and p_atom0[3] != 'CYS' and p_atom1[3] != 'MET' and \
                                                    p_atom1[3] != 'CYS':
                                                if 2.4 < length0 <= 3.2:
                                                    blijniy_4_2_A = [p_atom0, p_atom1, length0]

                                                if flag_d:
                                                    if 3.2 < length0 <= 6.4:
                                                        if blijniy_4_2_A == 0:  # to where check? (not in loop)
                                                            dalniy_8_0_A = [p_atom0, p_atom1, length0]

                                                if blijniy_4_2_A == 0 and p_atom0 not in alone_fg:
                                                    alone_fg.append(p_atom0)
                                                    alone = [p_atom0]  # not pair with blijniy
                                                    # print(f"Check repeat {alone}")

                                            if p_atom0[3] == 'MET' or p_atom0[3] == 'CYS':
                                                if p_atom1[3] != 'MET' and p_atom1[3] != 'CYS':
                                                    if 2.4 < length0 <= 3.5:
                                                        blijniy_4_2_A = [p_atom0, p_atom1, length0]

                                                    if flag_d:
                                                        if 3.5 < length0 <= 6.7:
                                                            if blijniy_4_2_A == 0:
                                                                dalniy_8_0_A = [p_atom0, p_atom1, length0]

                                                    if blijniy_4_2_A == 0 and p_atom0 not in alone_fg:
                                                        alone_fg.append(p_atom0)
                                                        alone = [p_atom0]
                                                        # print(f"Check repeat {alone}")

                                                if p_atom1[3] == 'MET' or p_atom1[3] == 'CYS':
                                                    if 2.4 < length0 <= 3.8:
                                                        blijniy_4_2_A = [p_atom0, p_atom1, length0, ]

                                                    if flag_d:
                                                        if 3.8 < length0 < 7.0:
                                                            if blijniy_4_2_A == 0:
                                                                dalniy_8_0_A = [p_atom0, p_atom1, length0]

                                                    if blijniy_4_2_A == 0 and p_atom0 not in alone_fg:
                                                        alone_fg.append(p_atom0)
                                                        alone = [p_atom0]
                                                        # print(f"Check repeat {alone}")

                                            if p_atom0[3] != 'MET' and p_atom0[3] != 'CYS':
                                                if p_atom1[3] == 'MET' or p_atom1[3] == 'CYS':
                                                    if 2.4 < length0 <= 3.5:
                                                        blijniy_4_2_A = [p_atom0, p_atom1, length0]

                                                    if flag_d:
                                                        if 3.5 < length0 <= 6.7:
                                                            if blijniy_4_2_A == 0:
                                                                dalniy_8_0_A = [p_atom0, p_atom1, length0]

                                                    if blijniy_4_2_A == 0 and p_atom0 not in alone_fg:
                                                        alone_fg.append(p_atom0)
                                                        alone = [p_atom0]
                                                        # print(f"Check repeat {alone}")

                                            if blijniy_4_2_A != 0 or dalniy_8_0_A != 0 or alone != 0:
                                                ion_pair = {'blijniy_4_2_A': blijniy_4_2_A,
                                                            'dalniy_8_0_A': dalniy_8_0_A,
                                                            'mix_alone': alone}
                                                ion_pairs.append(ion_pair)
                    if p_atom0 in p_atoms_copy:
                        p_atoms_copy.remove(p_atom0)
    return ion_pairs


# print('1_12_24')
