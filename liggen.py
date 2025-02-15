import datetime
import os
import sys
import json
import math
import subprocess
import src as lf


commands = ['q', 'h', 's', 'f', 'b', 'p', 'a', 'fi', 'fo', 'i', 'm', 'cpu']
commands_text = f"q (quit), h (help), s (start), b (bridges), m (mode), p (angle), a (atoms), fi (file input)," \
                f" cpu (load CPU), i (info)."
command_not = f"Command not found."
parameter_not = f"Parameter not found."
help_main = f"little description"


version = 0.1
cite = f"citation"
info_LigGen = f"______________________________________\n" \
              f"LigGen version: {version}.\n" \
              f"Developed by Aleksey Shustikov.\n" \
              f"Cite: {cite}.\n" \
              f"______________________________________\n"


def main_p(parrallelisation_file_full, process_time):
    print(f"\nSearching atom pairs... Please, wait at least two minutes.\nRemaining time is depend on amount of aminoacids in input file.\n")
    subprocess.run([sys.executable, parrallelisation_file_full])
    if os.path.exists('process_time'):
        os.rmdir(process_time)
    lf.data_end()


def liggen_cli():

    slash = lf.platform_get()
    cwd = os.getcwd()
    cfg_file_full = lf.configuration_make(file_name=True)

    processes_file_full = 0
    if os.path.exists(lf.joinprocesses_make(file_name=True)):
        lf.joinprocesses_make(remove=True)
    else:
        processes_file_full = lf.joinprocesses_make(file_name=True)
    parrallelisation_file_full = os.path.join(cwd, f"src{slash}" + 'for_parallel.py')

    if len(sys.argv) > 1:
        print(f"Only program name required.")

    if len(sys.argv) == 1:
        flag_hello = True
        while flag_hello:
            try:
                with open(cfg_file_full, 'r') as f:
                    f.read()
            except FileNotFoundError:
                lf.configuration_make()

            try:
                with open(processes_file_full, 'r') as f:
                    f.read()
            except FileNotFoundError:
                lf.joinprocesses_make()

            if not os.path.exists('Works'):
                os.mkdir('Works')
                os.mkdir('Tables')

            else:
                flag_hello = False
                flag_cfg = True
                while flag_cfg:
                    with open(cfg_file_full, 'r') as f:
                        contents = json.load(f)
                        if contents['file input'] is None:
                            flag_fi = True
                            while flag_fi:
                                fi_cfg = input(f"Enter absolute file path for .cif format file. Press 'q' to quit.\n")
                                fi_cfg = str(fi_cfg.lower())
                                if fi_cfg == 'q':
                                    lf.joinprocesses_make(remove=True)
                                    sys.exit()

                                file_path_0 = f"{fi_cfg}"
                                try:
                                    with open(file_path_0, "r") as fi:
                                        f_p_c = fi.read()
                                except FileNotFoundError:
                                    print(f"Absolute file path [{file_path_0}] not found.\n")
                                except OSError:
                                    pass
                                except EOFError:
                                    pass

                                else:
                                    with open(cfg_file_full, 'w') as f0:
                                        contents['file input'] = file_path_0
                                        json.dump(contents, f0)
                                        flag_fi = False

                        else:
                            flag_cfg = False
                            print(f"\n\tLiggen configuration file:")
                            for k, v in contents.items():
                                print(f"\tParameter [{k}]: {v}")

                            flag_first = True
                            while flag_first:
                                res_val = ['q', 's', 'c']
                                res0 = input(f"\nPress 's' to start a job or 'c' to change LigGen configuration file.\n"
                                             f"Press 'q' to quit.\n")
                                if res0 == 'q':
                                    lf.joinprocesses_make(remove=True)
                                    sys.exit()

                                if res0 == 's':
                                    main_p(parrallelisation_file_full=parrallelisation_file_full, process_time='process_time')
                                    flag_first = False

                                if res0 not in res_val:
                                    print(f"Value '{res0}' is not equal to 's' or 'c'.")

                                if res0 == 'c':
                                    flag_first = False
                                    flag_change = True
                                    while flag_change:
                                        command = input(f"\nEnter any command or press 'q' to exit.\n"
                                                        f"List of commands:\n{commands_text}\n")
                                        command = command.lower().strip()

                                        if command == "q":
                                            flag_change = False

                                        if command in commands:
                                            if command == "h":
                                                print(f"help page")

                                            if command == "s":
                                                main_p(parrallelisation_file_full=parrallelisation_file_full, process_time='process_time')
                                                flag_change = False

                                            if command == "b":
                                                flag_b = True
                                                while flag_b:
                                                    parameter_b = input(f"\nEnter any parameter for [{command}] command or 'q' to change command.\n")
                                                    if parameter_b == "q":
                                                        flag_b = False

                                                    elif parameter_b == "yes":
                                                        with open(cfg_file_full, 'w') as f0:
                                                            contents['bridges'] = parameter_b
                                                            json.dump(contents, f0)
                                                        flag_b = False

                                                    elif parameter_b == "no":
                                                        with open(cfg_file_full, 'w') as f0:
                                                            contents['bridges'] = parameter_b
                                                            json.dump(contents, f0)
                                                        flag_b = False

                                                    else:
                                                        print(f"Required only yes/no.")

                                            if command == "m":
                                                flag_m = True
                                                while flag_m:
                                                    parameter_b = input(f"\nEnter any parameter for [{command}] command or 'q' to change command.\n")
                                                    if parameter_b == "q":
                                                        flag_m = False

                                                    elif parameter_b == "prw":
                                                        with open(cfg_file_full, 'w') as f0:
                                                            contents['mode'] = parameter_b
                                                            json.dump(contents, f0)
                                                        flag_m = False

                                                    elif parameter_b == "prp":
                                                        with open(cfg_file_full, 'w') as f0:
                                                            contents['mode'] = parameter_b
                                                            json.dump(contents, f0)
                                                        flag_m = False

                                                    else:
                                                        print(f"Required only prw (predict relative to water) or prp "
                                                              f"(predict relative to protein).")

                                            if command == "p":
                                                flag_p = True
                                                while flag_p:
                                                    parameter_p = input(f"Enter an angle (2-90°) between neighbours points on a circle"
                                                                        f" (2 = 2°, 10 = 10°, etc) or press 'q' to change the command.\n")
                                                    if parameter_p == "q":
                                                        flag_p = False
                                                    else:
                                                        try:
                                                            parameter_p = math.fabs(int(parameter_p))
                                                        except ValueError:
                                                            print(f"Фальс! Parameter [{parameter_p}] must be a number.\n")
                                                        else:
                                                            if 2 <= parameter_p <= 90:
                                                                contents['precision'] = parameter_p
                                                                with open(cfg_file_full, 'w') as fp:
                                                                    json.dump(contents, fp)
                                                                flag_p = False
                                                            else:
                                                                print(f"{parameter_p} is not between 2 and 90.\n")

                                            if command == "cpu":
                                                flag_cpu = True
                                                while flag_cpu:
                                                    parameter_p = input(f"Enter desired нагрузку на процессор (load 30% by default) "
                                                                        f"or press 'q' to change the command.\n")
                                                    if parameter_p == "q":
                                                        flag_cpu = False
                                                    else:
                                                        try:
                                                            parameter_p = math.fabs(int(parameter_p))
                                                        except ValueError:
                                                            print(f"Фальс! Parameter [{parameter_p}] must be a number.")
                                                        else:
                                                            if 10 <= parameter_p <= 80:
                                                                parameter_p = f"{parameter_p}%"
                                                                contents['cpu count'] = parameter_p
                                                                with open(cfg_file_full, 'w') as fp:
                                                                    json.dump(contents, fp)
                                                                flag_cpu = False
                                                            else:
                                                                print(f"{parameter_p} is not between 10 and 80.\n")

                                            if command == "a":
                                                atoms = ['nitrogen', 'oxygen', 'sulphur']
                                                atoms_values = contents['atoms radii']

                                                flag_a = True
                                                for i, atom in enumerate(atoms):
                                                    while flag_a:
                                                        value = input(f"{atoms[i].title()} radii = ")
                                                        try:
                                                            value = math.fabs(float(value))
                                                        except ValueError:
                                                            print(f"Фальс! Parameter [{value}] must be a number.")
                                                        else:
                                                            if 2 < value < 5:
                                                                if i == 0:
                                                                    atoms_values['nitrogen'] = value

                                                                if i == 1:
                                                                    atoms_values['oxygen'] = value

                                                                if i == 2:
                                                                    atoms_values['sulphur'] = value

                                                                i += 1
                                                            else:
                                                                print(f"{atoms[i].title()} hydrate radii must be between 2-5 Angstrom.")
                                                            if i == len(atoms):
                                                                with open(cfg_file_full, 'w') as f1:
                                                                    json.dump(contents, f1)
                                                                flag_a = False

                                            if command == "fi":
                                                flag_fin = True
                                                print(f"\nCurrent file path to input:\n{contents['file input']}")
                                                while flag_fin:
                                                    parameter_f = input(f"\nType 'q' to change command or enter object name (.cif format).\n")
                                                    file_path = f"{parameter_f}"
                                                    if parameter_f == "q":
                                                        flag_fin = False
                                                    else:
                                                        try:
                                                            with open(file_path, "r") as fp:
                                                                contents_file = fp.read()
                                                        except FileNotFoundError:
                                                            print(f"\nAbsolute file path [{file_path}] not found.")
                                                        except OSError:
                                                            pass
                                                        else:
                                                            with open(cfg_file_full, 'w') as f0:
                                                                contents['file input'] = file_path
                                                                json.dump(contents, f0)
                                                            flag_fin = False

                                            if command == "fo":
                                                debug = False
                                                print(f"test, debug:{debug}.")
                                                if debug:
                                                    flag_fo = True
                                                    while flag_fo:
                                                        # while flag_fo: import pdb name, and date
                                                        parameter_fo = input(f"Enter 'd' if want to use a default file name output or 'c' "
                                                                             f"if want to use a custom file name output. Type q to change command.\n")
                                                        if parameter_fo == "q":
                                                            flag_fo = False
                                                        elif parameter_fo == 'd':
                                                            pdb_id = cp.get_pdb_id()
                                                            date_id = datetime.datetime.now()
                                                            print(pdb_id, date_id)
                                                            flag_fo = False
                                                        elif parameter_fo == 'c':
                                                            print(f"Used as a constant name for file.\n")
                                                            flag_fo = False
                                                        else:
                                                            print(f"Parameter '{parameter_fo}' is not equal to 'd' or 'c'.")

                                            if command == "i":
                                                flag_i = True
                                                while flag_i:
                                                    parameter_i = input(f"\nEnter 's' (session info) or 'l' (LigGen info). Type 'q' to change command.\n")
                                                    if parameter_i == "q":
                                                        flag_i = False
                                                    elif parameter_i == "s" or parameter_i == "l":
                                                        if parameter_i == "s":
                                                            print(f"\tPDB ID: {lf.fo_pdb_id()}")
                                                            for k, v in contents.items():
                                                                print(f"\tParameter [{k}]: {v}")
                                                            print()

                                                        if parameter_i == "l":
                                                            print(info_LigGen)
                                                    else:
                                                        print(f"Parameter [{parameter_i}] is not equal to s or l.")

                                        if command not in commands:
                                            print(f"{command_not}")

                                    # if flag_change != True:
                                    #     lf.joinprocesses_make(remove=True)


if __name__ == '__main__':
    liggen_cli()
