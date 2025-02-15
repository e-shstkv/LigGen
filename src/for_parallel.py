from multiprocessing.pool import Pool
import make_main as mm
import config as cfg


count_cpu = cfg.count_cpu()
lists_of_data = mm.process()


if __name__ == '__main__':
    with Pool(count_cpu) as pool:
        pool.map(mm.all_waters, lists_of_data)
