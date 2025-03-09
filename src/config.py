import os
import json
import platform


"""
File 'сonfig.py' takes info from CLI program and return prepared info to the 
script or program files (liggen_files [dir]):
"""


def platform_get():
    os_p = platform.system()
    os_slash = str
    if os_p == 'Windows':
        os_slash = f"\\"

    if os_p == 'Linux' or os_p == 'Darwin':
        os_slash = f"/"

    return os_slash


def configuration_make(print_result=None, file_name=None):
    """
    is file exist
    :return:Make 'configuration_liggen.JSON'
    """
    slash = platform_get()
    cwd = os.getcwd()
    cfg_file = 'configuration_liggen.JSON'
    cfg_file_full = os.path.join(cwd, f"src{slash}" + cfg_file)
    if print_result is None:
        try:
            with open(cfg_file_full, 'r') as f:
                cfg = json.load(f)
        except FileNotFoundError:
            print(f"\nMaking LigGen configuration file.")
            with open(cfg_file_full, 'w') as f:
                cfg = {'bridges': 'yes',
                       'mode': 'prp',
                       'precision': 10,
                       'atoms radii': {'nitrogen': 2.5, 'oxygen': 2.4, 'sulphur': 2.8},
                       'cpu count': '30.0%',
                       'file input': None}
                json.dump(cfg, f)
    if print_result:
        with open(cfg_file_full, 'r') as f:
            cfg = json.load(f)
            for k, v in cfg.items():
                print(f"{k}: {v}")
    if file_name:
        return cfg_file_full


def joinprocesses_make(file_name=None, remove=None):
    """
    is file exist
    :return:Make
    """
    slash = platform_get()
    cwd = os.getcwd()
    processes_file = 'joinprocesses.json'
    processes_file_full = os.path.join(cwd, f"src{slash}" + processes_file)

    if file_name is None and remove is None:
        try:
            with open(processes_file_full, 'r') as f:
                q = json.load(f)
        except FileNotFoundError:
            with open(processes_file_full, 'w') as f:
                content = {'oxygens': []}
                json.dump(content, f)
    if file_name:
        return processes_file_full

    if remove:
        os.remove(processes_file_full)


def configuration_read(key,):
    """
    Read LigGen configuration file.
    :return: value based on key
    """
    slash = platform_get()
    cwd = os.getcwd()
    cfg_file = 'configuration_liggen.JSON'
    cfg_file_full = os.path.join(cwd, f"src{slash}" + cfg_file)
    with open(cfg_file_full, 'r') as f:
        cfg = json.load(f)
        value = cfg[key]
        return value


def flag_distant():
    """
    :return:flag (take into account var(bridges -> True or None); flag == yes by default.
    """
    param = configuration_read('bridges')
    if param == 'yes':
        return True
    if param == 'no':
        return None


def count_cpu():
    """
    :return: % load of CPU
    """
    count = configuration_read('cpu count')
    cpu_n = []
    for q in count:
        if q == '.':
            break
        cpu_n.append(q)
    cpu_n = int(''.join(cpu_n))
    cpu_n = cpu_n / 100
    count_cpu = int(round((os.cpu_count() * cpu_n), 0))
    return count_cpu


def water_bond_parameter_radii(atom=None):
    """
    :param atom:
    :return: hydration radii of atom
    """
    atom_rad = configuration_read('atoms radii')
    if atom == 'O':
        return atom_rad['oxygen']
    if atom == 'N':
        return atom_rad['nitrogen']
    if atom == 'S':
        return atom_rad['sulphur']


def choose_precision():
    """
    Count of points on the circle based on angle between ATOM points
    :return: count of points on the circle
    """
    precision = int(round(360 / (configuration_read('precision')), 0))
    return precision


# print('14_01_25')
