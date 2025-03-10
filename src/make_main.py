import sympy
import math
from multiprocessing.pool import Pool
import os
import json
import time
import geometry_functions as gf
import content_parser as cp
import config as cfg


cwd = os.getcwd()
count_cpu = cfg.count_cpu()


def process():
    """
    make list of tasks
    :return: count of tasks, based on count CPU
    """
    ion_pairs = cp.find_ion_pairs_3_2()
    len_i_p = len(ion_pairs)
    copy_len_i_p = len_i_p
    count_div, slice_average, slice_last = 0, 0, 0
    flag_celoe = True
    while flag_celoe:
        div = copy_len_i_p / (count_cpu)
        if div % 1 == 0:
            slice_average = div
            slice_last = count_div
            flag_celoe = False

        if div % 1 != 0:
            copy_len_i_p += 1
            count_div += 1

    slice_list = [x for x in range(1 * count_cpu)]
    s_s = 0
    list_of_slices = []
    for q in slice_list:
        first = int(s_s)
        last = int(slice_average * q + slice_average)
        if q == slice_list[-1]:
            last = int(last - slice_last)

        do = ion_pairs[first: last]
        list_of_slices.append(do)
        s_s += slice_average

    return list_of_slices


def process_load(oxygen=None):
    """
    :param oxygen: one result (i)
    :return: results (i, i + 1, ...)
    """
    filename = cfg.joinprocesses_make(file_name=True)
    with open(filename, 'r+') as f:
        data_oxygen = json.load(f)
        data_oxygen['oxygens'].append(oxygen)
        f.seek(0)
        json.dump(data_oxygen, f)


def remaining_time(t_start=None, len_ion_pairs=None, n=None, last=None):
    end = round((time.time() - t_start) / 60, 1)
    time_k = len_ion_pairs / n
    remaining_time = round((end * time_k - end) * 2, 1)
    percent_done = round((n / len_ion_pairs) * 100, 1)
    promt = f"Time has passed: {end} minutes. Calculation is {percent_done}% complete.\nEstimated remaining time: {remaining_time} minutes.\n"
    if last is None:
        print(promt)
    if last:
        print(promt)
        print(f"Joining processess...\nWait at least 10 minutes.")


def process_time():
    if not os.path.exists('process_time'):
        os.mkdir('process_time')
        return 0
    else:
        return 1


def all_waters(ion_pairs=None):
    # t_process = time.time()
    qp = process_time()
    all_points = []
    len_ion_pairs = len(ion_pairs)
    if qp == 0:
        print(f"Found {len_ion_pairs * count_cpu} atom pairs. CPU count will be used: {count_cpu}\n")
    t_start = time.time()
    for n, i in enumerate(ion_pairs):
        if qp == 0:
            if n == int(len_ion_pairs * 0.01):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.05):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.1):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.2):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.4):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.6):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.8):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n)
            if n == int(len_ion_pairs * 0.99):
                remaining_time(t_start=t_start, len_ion_pairs=len_ion_pairs, n=n, last=True)

        for k, v in i.items():
            if k == 'blijniy_4_2_A' and v != 0:
                x0, y0, z0 = round(v[0][5][0], 1), round(v[0][5][1], 1), round(v[0][5][2], 1)
                x1, y1, z1 = round(v[1][5][0], 1), round(v[1][5][1], 1), round(v[1][5][2], 1)

                circle_fragments = gf.find_circle_fragments(atom0=v[0][1], x0=x0, y0=y0, z0=z0,
                                                            atom1=v[1][1], x1=x1, y1=y1, z1=z1, length=v[2])
                if circle_fragments:
                    for qc in circle_fragments:
                        all_points.append(qc)

            if k == 'dalniy_8_0_A' and v != 0:
                x0, y0, z0 = round(v[0][5][0], 1), round(v[0][5][1], 1), round(v[0][5][2], 1)
                x1, y1, z1 = round(v[1][5][0], 1), round(v[1][5][1], 1), round(v[1][5][2], 1)

                circle_fragments = gf.find_circle_fragments(atom0=v[0][1], x0=x0, y0=y0, z0=z0,
                                                            atom1=v[1][1], x1=x1, y1=y1, z1=z1, length=v[2])
                if circle_fragments:
                    for qc in circle_fragments:
                        all_points.append(qc)

    # t_process_end = time.time() - t_process
    # print(f"\nprocess end {t_process_end} seconds.")
    process_load(all_points)

# 09_03_25
