import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
from time import time
from multiprocessing import Process, Array
from math import ceil
import matplotlib.pyplot as plt



# Raise error if file not found

def check_file_path(file_path):
    if not file_path.is_file():
        raise FileNotFoundError(f'File {file_path} does not exist')



# Raise error if dir not found

def check_dir_path(dir_path):
    if not dir_path.is_dir():
        raise FileNotFoundError(f'Directory {dir_path} does not exist')



# Parse input file containing selected populations

def parse_selected_populations(file_path):
    with file_path.open(mode = 'r', encoding = 'utf-8') as file:
        content = file.read()
        pop_lines = [line.split() for line in content.splitlines() if not line.startswith('#')]
        sel_pops = [pop_line[0] for pop_line in pop_lines]

    return sel_pops



# Parse .ind input file containing population indices

def parse_populations(file_path):
    pops = defaultdict(list)
    num_pops = 0

    with file_path.open(mode = 'r', encoding = 'utf-8') as file:
        content = file.read()
        for index, line in enumerate(content.splitlines()):
            columns = line.split()
            pop_name = columns[-1]
            pops[pop_name].append(index)
            num_pops += 1

    return pops, num_pops



# Parse .snp input file containing allele names

def parse_snp_names(file_path):
    snp_names = []
    num_alleles = 0

    with file_path.open(mode = 'r', encoding = 'utf-8') as file:
        content = file.read()
        for line in content.splitlines():
            columns = line.split()
            snp_names.append(columns[0])
            num_alleles += 1

    return snp_names, num_alleles



# Count number of rows and columns in .geno input file

def geno_table_shape(file_path):
    num_rows = 0
    num_columns = []

    with file_path.open(mode = 'r', encoding = 'utf-8') as file:
        for row in file:
            num_rows += 1
            num_columns.append(len(row) - 1)

    if not all(nc == num_columns[0] for nc in num_columns):
        raise Exception(f'Error in {file_path}: Not all rows are of equal number of columns')

    return num_rows, num_columns[0]



# Check all input files consistency

def check_input_files(num_alleles, num_rows, num_populations, num_columns):
    if num_alleles != num_rows:
        raise Exception(f'Number of alleles ({num_alleles}) in .snp file is not equal to number of rows ({num_rows}) in .geno file')
    if num_populations != num_columns:
        raise Exception(f'Number of populations ({num_populations}) in .ind file is not equal to number of columns ({num_columns}) in .geno file')



# Compute frequencies given list of alleles

def allele_frequency(alleles):
    freq = 0
    num_alleles = 0

    for a in alleles:
        if a != 9:
            freq += (2 - a) / 2
            num_alleles += 1

    if num_alleles == 0:
        return -1

    return freq / num_alleles



# Compute frequencies of a population

def population_allele_frequencies(file_path, pop_indices, allele_freqs):
    with file_path.open(mode = 'r', encoding = 'utf-8') as file:
        for index, row in enumerate(file):
            allele_freqs[index] = allele_frequency([int(row[i]) for i in pop_indices])



# Parallel compute frequencies of all populations

def parallel_compute_populations_frequencies(num_procs, num_alleles, file_path, pop_indices, sel_pops):
    num_computations = len(pop_indices)
    batch_size = ceil(num_computations / num_procs)
    index = 0

    print(f'Computing {num_alleles} frequencies per population for {num_computations} populations in {batch_size} batches of {num_procs} parallel processes...')

    allele_freqs = [Array('d', num_alleles) for i in range(num_computations)]

    for nb in range(batch_size):
        procs = []
        computing_pops = []

        for nt in range(num_procs):
            if index < num_computations:
                p = Process(target = population_allele_frequencies, args = (file_path, pop_indices[index], allele_freqs[index]))
                procs.append(p)
                p.start()
                computing_pops.append(sel_pops[index])
            index += 1

        print("Computing populations:", ' '.join(computing_pops))

        for p in procs:
            p.join()

    return [np.array(freqs) for freqs in allele_freqs]



# Select allele indices for which all frequencies have a values of 9

def invalid_allele_indices(allele_freqs):
    invalid_indices = []

    for freqs in allele_freqs:
        for index, freq in enumerate(freqs):
            if freq == -1:
                invalid_indices.append(index)

    invalid_indices = np.unique(np.array(invalid_indices, dtype = int))

    return invalid_indices



# Select only valid indices

def remove_invalid_alleles(allele_freqs, indices):
    for i in range(len(allele_freqs)):
        allele_freqs[i] = np.delete(allele_freqs[i], indices)



# Computation of alpha pre JL

def mixing_coefficient_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs):
    parent_diff = parent1_freqs - parent2_freqs
    return np.dot(hybrid_freqs - parent2_freqs, parent_diff) / np.dot(parent_diff, parent_diff)



# Computation of admixture angle pre JL

def admixture_angle_pre_jl(hybrid_freqs, parent1_freqs, parent2_freqs):
    xa = hybrid_freqs - parent1_freqs
    xb = hybrid_freqs - parent2_freqs

    cosine = np.dot(xa, xb) / np.sqrt(np.dot(xa, xa) * np.dot(xb, xb))
    angle = np.arccos(cosine)
    percentage = angle / np.pi

    return cosine, angle * 180 / np.pi, percentage



# Computation of f3

def f3(hybrid_freqs, parent1_freqs, parent2_freqs):
    num_alleles = hybrid_freqs.size
    return np.dot(hybrid_freqs - parent1_freqs, hybrid_freqs - parent2_freqs) / num_alleles



# Computation of f4 prime

def f4_prime(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):
    num_aux_pops = len(aux_freqs)
    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)

    f4ab_prime = np.zeros(num_pairs)
    f4xb_prime = np.zeros(num_pairs)

    ab = parent1_freqs - parent2_freqs
    xb = hybrid_freqs - parent2_freqs

    index = 0

    for i in range(num_aux_pops):
        for j in range(i + 1, num_aux_pops):
            ij = aux_freqs[i] - aux_freqs[j]
            norm_ij = np.linalg.norm(ij)
            f4ab_prime[index] = np.dot(ab, ij) / norm_ij
            f4xb_prime[index] = np.dot(xb, ij) / norm_ij
            index += 1

    return f4ab_prime, f4xb_prime



# Computation of f4 standard

def f4_std(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):
    num_aux_pops = len(aux_freqs)
    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)

    f4ab_std = np.zeros(num_pairs)
    f4xb_std = np.zeros(num_pairs)

    ab = parent1_freqs - parent2_freqs
    xb = hybrid_freqs - parent2_freqs

    index = 0

    for i in range(num_aux_pops):
        for j in range(i + 1, num_aux_pops):
            ij = aux_freqs[i] - aux_freqs[j]
            f4ab_std[index] = np.dot(ab, ij)
            f4xb_std[index] = np.dot(xb, ij)
            index += 1

    num_alleles = hybrid_freqs.size

    return f4ab_std / num_alleles, f4xb_std / num_alleles



# Least squares fit

def least_squares(x, y):
    dim = len(x)

    A = np.vstack([x, np.zeros(dim)]).T
    alpha = np.linalg.lstsq(A, y)[0][0]

    Q = 0
    for i in range(dim):
        Q += (y[i] - alpha * x[i]) ** 2

    x_avg = 0
    for i in range(dim):
        x_avg += x[i]
    x_avg /= dim

    x_dev = 0
    for i in range(dim):
        x_dev += (x[i] - x_avg) ** 2

    s_alpha = np.sqrt(Q / ((dim - 2) * x_dev))
    t = 1.98

    error = s_alpha * t

    return alpha, error



# Computation of admixture angle post JL

def admixture_angle_post_jl(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):
    num_aux_pops = len(aux_freqs)
    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)

    xa = hybrid_freqs - parent1_freqs
    xb = hybrid_freqs - parent2_freqs

    sum1 = 0
    sum2 = 0
    sum3 = 0


    for i in range(num_aux_pops):
        for j in range(i + 1, num_aux_pops):
            ij = aux_freqs[i] - aux_freqs[j]

            xaij = np.dot(xa, ij)
            xbij = np.dot(xb, ij)
            ijij = np.dot(ij, ij)

            sum1 += xaij * xbij / ijij
            sum2 += (xaij ** 2) / ijij
            sum3 += (xbij ** 2) / ijij

    cosine = sum1 / np.sqrt(sum2 * sum3)
    angle = np.arccos(cosine)
    percentage = angle / np.pi

    return cosine, angle * 180 / np.pi, percentage



# Computation of f4-ratio

def f4_ratio(hybrid_freqs, parent1_freqs, parent2_freqs, aux_freqs):
    num_aux_pops = len(aux_freqs)
    num_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)

    xb = hybrid_freqs - parent2_freqs
    ab = parent1_freqs - parent2_freqs

    alpha = np.zeros(num_pairs)

    index = 0

    for i in range(num_aux_pops):
        for j in range(i + 1, num_aux_pops):
            ij = aux_freqs[i] - aux_freqs[j]
            alpha[index] = np.dot(xb, ij) / np.dot(ab, ij)
            index += 1

    alpha_01 = alpha[(alpha >= 0) & (alpha <= 1)]
    alpha_avg = np.average(alpha_01)
    alpha_std_dev = np.std(alpha_01) * 1.96
    alpha_hist = np.histogram(alpha, 20)

    return alpha, alpha_avg, alpha_std_dev, alpha_hist, alpha_01.size



# Ploat a fit

def plot_fit(x, y, alpha, title, xlabel, ylabel):
    fig, ax = plt.subplots()

    ax.set_title(title)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.plot(x, y, '.')
    ax.plot(x, alpha * x)

    plt.show()



# Plot a histogram

def plot_histogram(histogram, title, xlabel, ylabel):
    counts = histogram[0]
    edges = histogram[1]

    fig, ax = plt.subplots()

    ax.set_title(title)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.bar(edges[:-1], counts, width = np.diff(edges), edgecolor = 'black', align = 'edge')

    plt.show()



# Save frequencies

def save_population_allele_frequencies(out_dir_path, populations, allele_freqs):
    file_path = out_dir_path / 'frequencies.dat'

    with file_path.open(mode = 'w', encoding = 'utf-8') as file:
        pops_width = max([len(name) for name in populations])
        prec = 6
        col_width = max(prec + 7, pops_width)

        headers_format = ' '.join([f'{{{i}:^{col_width}}}' for i, pop in enumerate(populations)])
        headers = headers_format.format(*populations)
        file.write(headers + '\n')

        row_format = ' '.join([f'{{{i}: {col_width}.{prec}E}}' for i, pop in enumerate(populations)])

        num_alleles = allele_freqs[0].size

        for allele_index in range(num_alleles):
            row = [freqs[allele_index] for freqs in allele_freqs]
            file.write(row_format.format(*row) + '\n')



# Save f4 points

def save_f4_points(out_dir_path, f4ab_prime, f4xb_prime, f4ab_std, f4xb_std, alpha, sel_pops):
    file_path = out_dir_path / 'f4-points.dat'

    with file_path.open(mode = 'w', encoding = 'utf-8') as file:
        aux_pops_width = max([len(name) for name in sel_pops[3:]])
        prec = 6
        col_width = prec + 7

        headers = '{0:^{col_width}} {1:^{col_width}} {2:^{col_width}} {3:^{col_width}} {4:^{col_width}} {5:^{aux_pops_width}} {6:^{aux_pops_width}}'.format('f4primeAB', 'f4primeXB', 'f4AB', 'f4XB', 'f4-ratio', 'Aux1', 'Aux2', col_width = col_width, aux_pops_width = aux_pops_width)
        file.write(headers + '\n')

        num_aux_pops = len(sel_pops[3:])

        index = 0

        for i in range(num_aux_pops):
            for j in range(i + 1, num_aux_pops):
                row = '{0: {col_width}.{prec}E} {1: {col_width}.{prec}E} {2: {col_width}.{prec}E} {3: {col_width}.{prec}E} {4: {col_width}.{prec}E} {5:{aux_pops_width}} {6:{aux_pops_width}}'.format(f4ab_prime[index], f4xb_prime[index], f4ab_std[index], f4xb_std[index], alpha[index], sel_pops[3:][i], sel_pops[3:][j], prec = prec, col_width = col_width, aux_pops_width = aux_pops_width)
                file.write(row + '\n')
                index += 1



# Save results

def save_admixture_data(out_dir_path, num_alleles, num_valid_alleles, sel_pops, cosine_pre_jl, angle_pre_jl, percentage_pre_jl, cosine_post_jl, angle_post_jl, percentage_post_jl, alpha_pre_jl, alpha, alpha_error, alpha_std, alpha_std_error, alpha_ratio_avg, alpha_ratio_std_dev, num_cases, f3_test):
    file_path = out_dir_path / 'results.dat'

    with file_path.open(mode = 'w', encoding = 'utf-8') as file:
        num_aux_pops = len(sel_pops[3:])
        num_aux_pairs = int(num_aux_pops * (num_aux_pops - 1) / 2)

        prec = 6

        file.write(f'Admixture model: {sel_pops[0]} = {sel_pops[1]} + {sel_pops[2]}\n')
        file.write(f'SNPs employed: {num_valid_alleles} / {num_alleles}\n')
        file.write(f'Auxiliary populations: {num_aux_pops}\n')
        file.write(f'Auxiliary pairs: {num_aux_pairs}\n')
        file.write(f'Cos pre-JL:  {cosine_pre_jl:7.4f} ---> Angle pre-JL:  {angle_pre_jl:7.2f} deg vs 180 deg: {percentage_pre_jl:.1%}\n')
        file.write(f'Cos post-JL: {cosine_post_jl:7.4f} ---> Angle post-JL: {angle_post_jl:7.2f} deg vs 180 deg: {percentage_post_jl:.1%}\n')
        file.write(f'Alpha pre-JL:     {alpha_pre_jl:6.4f}\n')
        file.write(f'Alpha post-JL:    {alpha:6.4f} +/- {alpha_error:6.4f} (95% CI) (f4-prime, renormalized)\n')
        file.write(f'Alpha NR post-JL: {alpha_std:6.4f} +/- {alpha_std_error:6.4f} (95% CI) (f4, standard)\n')
        file.write(f'f4-ratio average if [0, 1]: {alpha_ratio_avg:6.4f} +/- {alpha_ratio_std_dev:6.4f} (95% CI), {num_cases} cases\n')
        file.write(f'Standard admixture test: f3(c1, c2; x) < 0 ? {f3_test:8.6f}')

def save_timings(out_dir_path, exec_times):
    file_path = out_dir_path / 'timings.dat'

    with file_path.open(mode = 'w', encoding = 'utf-8') as file:
        for key, value in exec_times.items():
            file.write(f'{key}: {value}\n')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Mixtum: The geometry of admixture in population genetics')
    parser.add_argument('--geno', type = str, required = True, help = 'Path of .geno file')
    parser.add_argument('--ind', type = str, required = True, help = 'Path of .ind file')
    parser.add_argument('--snp', type = str, required = True, help = 'Path of .snp file')
    parser.add_argument('--pops', type = str, required = True, help = 'Path of selected populations file')
    parser.add_argument('--outdir', type = str, required = True, help = 'Path of output dir')
    parser.add_argument('--nprocs', type = int, default = 1, help = 'Number of parallel computation processes')

    args = parser.parse_args()

    num_procs = args.nprocs

    geno_file_path = Path(args.geno)
    ind_file_path = Path(args.ind)
    snp_file_path = Path(args.snp)
    pops_file_path = Path(args.pops)
    out_dir_path = Path(args.outdir)

    check_file_path(geno_file_path)
    check_file_path(ind_file_path)
    check_file_path(snp_file_path)
    check_file_path(pops_file_path)
    check_dir_path(out_dir_path)

    exec_times = dict()
    exec_times['num_processes'] = num_procs

    t1 = time()

    sel_pops = parse_selected_populations(pops_file_path)
    populations, num_populations = parse_populations(ind_file_path)
    snp_names, num_alleles = parse_snp_names(snp_file_path)

    t2 = time()
    exec_times['parse_inputs'] = t2 - t1

    t1 = time()

    num_geno_rows, num_geno_columns = geno_table_shape(geno_file_path)
    check_input_files(num_alleles, num_geno_rows, num_populations, num_geno_columns)

    t2 = time()
    exec_times['check_inputs'] = t2 - t1

    t1 = time()

    pop_indices = [populations[pop] for pop in sel_pops]
    allele_freqs = parallel_compute_populations_frequencies(num_procs, num_alleles, geno_file_path, pop_indices, sel_pops)

    t2 = time()
    exec_times['compute_frequencies'] = t2 - t1

    t1 = time()

    invalid_indices = invalid_allele_indices(allele_freqs)
    num_valid_alleles = num_alleles - invalid_indices.size
    remove_invalid_alleles(allele_freqs, invalid_indices)

    t2 = time()
    exec_times['purge_indices'] = t2 - t1

    t1 = time()

    alpha_pre_jl = mixing_coefficient_pre_jl(allele_freqs[0], allele_freqs[1], allele_freqs[2])

    cosine_pre_jl, angle_pre_jl, percentage_pre_jl = admixture_angle_pre_jl(allele_freqs[0], allele_freqs[1], allele_freqs[2])

    f3_test = f3(allele_freqs[0], allele_freqs[1], allele_freqs[2])

    f4ab_prime, f4xb_prime = f4_prime(allele_freqs[0], allele_freqs[1], allele_freqs[2], allele_freqs[3:])
    alpha, alpha_error = least_squares(f4ab_prime, f4xb_prime)

    f4ab_std, f4xb_std = f4_std(allele_freqs[0], allele_freqs[1], allele_freqs[2], allele_freqs[3:])
    alpha_std, alpha_std_error = least_squares(f4ab_std, f4xb_std)

    cosine_post_jl, angle_post_jl, percentage_post_jl = admixture_angle_post_jl(allele_freqs[0], allele_freqs[1], allele_freqs[2], allele_freqs[3:])

    alpha_ratio, alpha_ratio_avg, alpha_ratio_std_dev, alpha_ratio_hist, num_cases = f4_ratio(allele_freqs[0], allele_freqs[1], allele_freqs[2], allele_freqs[3:])

    t2 = time()
    exec_times['compute_results'] = t2 - t1

    plot_fit(f4ab_prime, f4xb_prime, alpha, f'Renormalized admixture: {sel_pops[0]} = alpha {sel_pops[1]} + (1 - alpha) {sel_pops[2]}', f"f4'({sel_pops[1]}, {sel_pops[2]}; i, j)", f"f4'({sel_pops[0]}, {sel_pops[2]}; i, j)")
    plot_fit(f4ab_std, f4xb_std, alpha_std, f'Standard admixture: {sel_pops[0]} = alpha {sel_pops[1]} + (1 - alpha) {sel_pops[2]}', f"f4({sel_pops[1]}, {sel_pops[2]}; i, j)", f"f4({sel_pops[0]}, {sel_pops[2]}; i, j)")

    plot_histogram(alpha_ratio_hist, f'{sel_pops[0]} = alpha {sel_pops[1]} + (1 - alpha) {sel_pops[2]}', 'f4 ratio', 'Counts')

    t1 = time()

    save_population_allele_frequencies(out_dir_path, sel_pops, allele_freqs)
    save_f4_points(out_dir_path, f4ab_prime, f4xb_prime, f4ab_std, f4xb_std, alpha_ratio, sel_pops)
    save_admixture_data(out_dir_path, num_alleles, num_valid_alleles, sel_pops, cosine_pre_jl, angle_pre_jl, percentage_pre_jl, cosine_post_jl, angle_post_jl, percentage_post_jl, alpha_pre_jl, alpha, alpha_error, alpha_std, alpha_std_error, alpha_ratio_avg, alpha_ratio_std_dev, num_cases, f3_test)

    t2 = time()
    exec_times['save_outputs'] = t2 - t1

    save_timings(out_dir_path, exec_times)
