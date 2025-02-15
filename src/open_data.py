import os
import math
import json
from content_parser import get_pdb_id
import config as cfg


def fo_pdb_id():
    pdb_id = get_pdb_id()
    return pdb_id


def open_processes(joinprocesses):
    oxygens = []
    with open(joinprocesses, 'r') as f_n:
        content = json.load(f_n)
        for kline, vline in content.items():
            if kline == 'oxygens' and len(vline) != 0:
                for water0 in vline:
                    for water1 in water0:

                        sort_value = water1[1] - 1
                        coord = water1[0]
                        oxygen = [sort_value, coord]
                        oxygens.append(oxygen)
    return oxygens


def find_relevant_oxygens(oxygens=None):
    """
    Works with sorted oxygens by p-value
    :returns relevant (by p-value) oxygens;
    """
    oxygens_copy = oxygens[:]
    for i0, oxygen0 in enumerate(oxygens):
        coord0 = oxygen0[1]
        x0, y0, z0 = coord0[0], coord0[1], coord0[2]
        for i1, oxygen1 in enumerate(oxygens_copy):
            if oxygen0 != oxygen1:
                coord1 = oxygen1[1]
                x1, y1, z1 = coord1[0], coord1[1], coord1[2]
                length0 = math.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2 + (z0 - z1) ** 2)
                if length0 < 2.4:
                    if oxygen1 in oxygens:
                        oxygens.remove(oxygen1)
    return oxygens


def filename_check(file_check=None, pdb_id=None, dir=None, type=None):
    """
    :return: None or new filename
    """
    name = dir[0:-1]
    works_dir = os.listdir(dir)
    if pdb_id:
        pdb_works = []
        for f_n in works_dir:
            f_n_s = f_n.split()
            if len(f_n_s) >= 5 and f_n_s[2] == pdb_id:
                pdb_works.append(int(f_n_s[3]))

        if len(pdb_works) != 0:
            pdb_works.sort(reverse=True)
            update = int(pdb_works[0] + 1)
            update_name = f'{name}_predicted for {pdb_id} {update} times{type}'
            return update_name

        else:
            return file_check

    else:
        if file_check not in works_dir:
            print(f"{file_check} not in works dir")
            return file_check

        else:
            pdb_works = []
            for f_n in works_dir:
                f_n_s = f_n.split()
                if len(f_n_s) >= 5 and f_n_s[2] == str(pdb_id):
                    pdb_works.append(int(f_n_s[3]))

            if len(pdb_works) != 0:
                pdb_works.sort(reverse=True)
                update = int(pdb_works[0] + 1)
                update_name = f'{name}_predicted for {pdb_id} {update} times{type}'
                return update_name


def data_end():
    print(f"\nSorting predicted points...")
    joinprocesses = cfg.joinprocesses_make(file_name=True)
    oxygens = open_processes(joinprocesses=joinprocesses)
    oxygens.sort(key=lambda coordinate: coordinate[0], reverse=True)
    pdb_id = fo_pdb_id()
    cwd = os.getcwd()

    work_file = f'Work_predicted for {pdb_id} 1 times.cif'
    table_file = f'Table_predicted for {pdb_id} 1 times.txt'

    work_file_check = filename_check(file_check=work_file, pdb_id=pdb_id, dir='Works', type='.cif')
    fullname_w = os.path.join(cwd, "Works\\" + work_file_check)
    table_file_check = filename_check(file_check=table_file, pdb_id=pdb_id, dir='Tables', type='.txt')
    fullname_t = os.path.join(cwd, "Tables\\" + table_file_check)

    with open(fullname_w, 'w') as f:
        reduce_oxygens = find_relevant_oxygens(oxygens=oxygens)
        for i, oxygen in enumerate(reduce_oxygens, start=1):
            water_percent = round(float(oxygen[0]), 3)
            water_xyz = oxygen[1]
            x, y, z = round(water_xyz[0], 3), round(water_xyz[1], 3), round(water_xyz[2], 3)
            line = f"HETATM {water_percent} O  O   . HOH H 5 .   ? {x}  {y}  {z} 1.00 00.00  ? {i} HOH 0 O   1 "
            f.write(line + '\n')

    with open(fullname_t, 'w') as f:
        reduce_oxygens = find_relevant_oxygens(oxygens=oxygens)
        shapka = f"№    | score "
        f.write(shapka + '\n')
        for i, oxygen in enumerate(reduce_oxygens, start=1):
            water_percent = str(round(float(oxygen[0]), 3))
            i0 = str(i).ljust(5)
            wp_param = 5
            if water_percent[0] != '-':
                wp_param = 6
            wp = water_percent.ljust(wp_param)
            line = f"{i0}| {wp}"
            f.write(line + '\n')
    os.remove(joinprocesses)

    print(f"\nSorting has been finished.\nSee results in files:\n{fullname_w}\n{fullname_t}")



