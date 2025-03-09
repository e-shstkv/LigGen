import time
import math
import sympy
import content_parser as cp
import config as cfg


def find_circle_fragments(atom0=None, x0=None, y0=None, z0=None, atom1=None, x1=None, y1=None, z1=None, length=None):
    """
    find points on the circle, formed by intersection of two spheres;
    points (i, i + 1 ...) separated by constant length;
    :return: list of points [based on precision].
    """
    # time_start = time.time()

    # Sympy boost ->
    coords_0 = [x0, y0, z0]
    coords_1 = [x1, y1, z1]
    x0e, y0e, z0e, x1e, y1e, z1e = x0, y0, z0, x1, y1, z1
    div_x, div_y, div_z = 0, 0, 0
    for i, point_0, point_1 in zip(enumerate(coords_0), coords_0, coords_1):
        if point_0 >= 10 or point_1 >= 10:
            if point_0 >= point_1:
                if i[0] == 0:
                    x0e, x1e = x0 - point_0, x1 - point_0
                    div_x = point_0

                if i[0] == 1:
                    y0e, y1e = y0 - point_0, y1 - point_0
                    div_y = point_0

                if i[0] == 2:
                    z0e, z1e = z0 - point_0, z1 - point_0
                    div_z = point_0

            if point_0 < point_1:
                if i[0] == 0:
                    x0e, x1e = x0 - point_1, x1 - point_1
                    div_x = point_1

                if i[0] == 1:
                    y0e, y1e = y0 - point_1, y1 - point_1
                    div_y = point_1

                if i[0] == 2:
                    z0e, z1e = z0 - point_1, z1 - point_1
                    div_z = point_1

        if point_0 <= -10 or point_1 <= -10:
            if point_0 >= point_1:
                if i[0] == 0:
                    x0e, x1e = x0 - point_1, x1 - point_1
                    div_x = point_1

                if i[0] == 1:
                    y0e, y1e = y0 - point_1, y1 - point_1
                    div_y = point_1

                if i[0] == 2:
                    z0e, z1e = z0 - point_1, z1 - point_1
                    div_z = point_1

            if point_0 < point_1:
                if i[0] == 0:
                    x0e, x1e = x0 - point_0, x1 - point_0
                    div_x = point_0

                if i[0] == 1:
                    y0e, y1e = y0 - point_0, y1 - point_0
                    div_y = point_0

                if i[0] == 2:
                    z0e, z1e = z0 - point_0, z1 - point_0
                    div_z = point_0
    # Sympy boost <-

    p_atoms = cp.protein_atoms()
    polar_atoms = ['N', 'O', 'S']
    parametr_sort = True
    precision = cfg.choose_precision()
    atom00 = cfg.water_bond_parameter_radii(atom=atom0)
    atom11 = cfg.water_bond_parameter_radii(atom=atom1)
    if atom0:
        if length / (atom00 + atom11) < 0.95:  # have solutions   && (5 % RMSD)
            x2, y2, z2 = sympy.symbols('x y z', real=True)
            eq2 = sympy.Eq((x2 - x0e) ** 2 + (y2 - y0e) ** 2, atom00 ** 2 - (z2 - z0e) ** 2)
            eq3 = sympy.Eq((x2 - x1e) ** 2 + (y2 - y1e) ** 2, atom11 ** 2 - (z2 - z1e) ** 2)
            eq4 = sympy.Eq(z2, (z0e + z1e) / 2)
            answer1 = sympy.solve([eq2, eq3, eq4])
            answer0 = []

            for q in answer1[1:2]:
                for k, i in q.items():
                    if k == x2:
                        answer0.append(i)

                    if k == y2:
                        answer0.append(i)

                    if k == z2:
                        answer0.append(i)

            V1 = sympy.Matrix([x1, y1, z1])
            E = sympy.Matrix([(x0 - x1) / length, (y0 - y1) / length, (z0 - z1) / length])
            if len(answer0) > 1:
                A = sympy.Matrix([(answer0[0] + div_x) - x1, (answer0[1] + div_y) - y1, (answer0[2] + div_z) - z1])
                x, y, z = E[0], E[1], E[2]
                coordinates = []
                i = 1
                while i <= precision:
                    t = (2 * math.pi / precision) * i
                    M = sympy.Matrix(
                        [[(math.cos(t) + (1 - math.cos(t)) * x ** 2), ((1 - math.cos(t)) * x * y - math.sin(t) * z),
                          ((1 - math.cos(t)) * x * z + math.sin(t) * y)],
                         [((1 - math.cos(t)) * x * y + math.sin(t) * z), (math.cos(t) + (1 - math.cos(t)) * y ** 2),
                          ((1 - math.cos(t)) * y * z - math.sin(t) * x)],
                         [((1 - math.cos(t)) * x * z - math.sin(t) * y), ((1 - math.cos(t)) * y * z + math.sin(t) * x),
                          (math.cos(t) + (1 - math.cos(t)) * z ** 2)]])
                    coord_sympy = M * A + V1
                    xa, ya, za = float(coord_sympy[0]), float(coord_sympy[1]), float(coord_sympy[2])
                    coordinate = [xa, ya, za]
                    coordinates.append(coordinate)
                    i += 1

                exclude = []
                coordinates_and_values = []

                for i0, coordinate0 in enumerate(coordinates):
                    x3, y3, z3 = coordinate0[0], coordinate0[1], coordinate0[2]
                    for p_atom0 in p_atoms:
                        x4, y4, z4 = round(p_atom0[5][0], 1), round(p_atom0[5][1], 1), round(p_atom0[5][2], 1)
                        length0 = math.sqrt((x3 - x4) ** 2 + (y3 - y4) ** 2 + (z3 - z4) ** 2)
                        if length0 < 2.4:
                            exclude.append(coordinate0)

                for coordinate0 in coordinates:
                    if coordinate0 not in exclude:
                        polar, carbon, f_polar, f_carbon = 1, 1, 1, 1
                        x3, y3, z3 = coordinate0[0], coordinate0[1], coordinate0[2]
                        for p_atom0 in p_atoms:
                            x4, y4, z4 = round(p_atom0[5][0], 1), round(p_atom0[5][1], 1), round(p_atom0[5][2], 1)
                            length0 = math.sqrt((x3 - x4) ** 2 + (y3 - y4) ** 2 + (z3 - z4) ** 2)

                            if 2.4 <= length0 <= 3.2:
                                if p_atom0[1] == 'C':
                                    carbon += 1

                                if p_atom0[1] in polar_atoms:
                                    polar += 1

                            if 3.2 < length0 <= 5.0:
                                if p_atom0[1] == 'C':
                                    d0 = length0 - 3.2
                                    d1 = d0 / 1.8
                                    d2 = 1 - d1
                                    f_carbon += d2

                                if p_atom0[1] in polar_atoms:
                                    d0 = length0 - 3.2
                                    d1 = d0 / 1.8
                                    d2 = 1 - d1
                                    f_polar += d2

                        if parametr_sort:
                            int_val = polar / carbon
                            float_val = f_polar / f_carbon
                            val_env = int_val + float_val
                            pnt = [coordinate0, val_env]
                            coordinates_and_values.append(pnt)

                # end = (time.time() - time_start)
                # print(f"seconds per circle: {end}")

                if len(coordinates_and_values) != 0:
                    return coordinates_and_values

                else:
                    return None

        if length / (atom00 + atom11) > 0.95:  # has a solution
            exclude = []
            avrg_coord = (x1 + x0) / 2, (y1 + y0) / 2, (z0 + z1) / 2
            x5, y5, z5 = avrg_coord[0], avrg_coord[1], avrg_coord[2]
            polar, carbon, f_polar, f_carbon = 1, 1, 1, 1
            for p_atom0 in p_atoms:
                x6, y6, z6 = round(p_atom0[5][0], 1), round(p_atom0[5][1], 1), round(p_atom0[5][2], 1)
                length0 = math.sqrt((x5 - x6) ** 2 + (y5 - y6) ** 2 + (z5 - z6) ** 2)
                if length0 < 2.4:
                    exclude.append(avrg_coord)

            if avrg_coord not in exclude:
                for p_atom0 in p_atoms:
                    x6, y6, z6 = round(p_atom0[5][0], 1), round(p_atom0[5][1], 1), round(p_atom0[5][2], 1)
                    length0 = math.sqrt((x5 - x6) ** 2 + (y5 - y6) ** 2 + (z5 - z6) ** 2)

                    if 2.4 <= length0 <= 3.2:
                        if p_atom0[1] == 'C':
                            carbon += 1

                        if p_atom0[1] in polar_atoms:
                            polar += 1

                    if 3.2 < length0 <= 5.0:
                        if p_atom0[1] == 'C':
                            d0 = length0 - 3.2
                            d1 = d0 / 1.8
                            d2 = 1 - d1
                            f_carbon += d2

                        if p_atom0[1] in polar_atoms:
                            d0 = length0 - 3.2
                            d1 = d0 / 1.8
                            d2 = 1 - d1
                            f_polar += d2

            if avrg_coord in exclude:
                return None

            if parametr_sort:
                int_val = polar / carbon
                float_val = f_polar / f_carbon
                val_env = int_val + float_val
                pnt = [avrg_coord, val_env]
                return [pnt]

            else:
                return None


# print('02_01_25')
