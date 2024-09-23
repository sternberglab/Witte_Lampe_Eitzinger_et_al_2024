import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import gzip

def read_bedgraph(bedgraph_file):
    positions = {}
    coverage_values = {}
    
    with gzip.open(bedgraph_file, 'rt') as f:
        for line in f:
            if line.startswith('track') or line.startswith('#'):
                continue
            chrom, start, end, coverage = line.strip().split()
            start, end, coverage = int(start), int(end), float(coverage)
            
            if chrom not in positions:
                positions[chrom] = []
                coverage_values[chrom] = []
            
            positions[chrom].extend(range(start, end))
            coverage_values[chrom].extend([coverage] * (end - start))
    
    return positions, coverage_values

def bin_data(positions, coverage_values, bin_size, mask_size=0):
    binned_positions = {}
    binned_coverages = {}
    
    for chrom in positions.keys():
        if chrom == 'chrY':
            continue  # Skip chromosome Y
        
        binned_positions[chrom] = []
        binned_coverages[chrom] = []
        
        if not positions[chrom] or not coverage_values[chrom]:
            print(f"{chrom}: No data available.")
            continue
        
        current_bin_start = positions[chrom][0]
        current_bin_end = current_bin_start + bin_size
        bin_sum = 0
        bin_count = 0

        for pos, cov in zip(positions[chrom], coverage_values[chrom]):
            if chrom == 'chr1' and pos < mask_size:  # Mask the first `mask_size` bases of chromosome 1
                continue
            
            if pos < current_bin_end:
                bin_sum += cov
                bin_count += 1
            else:
                if bin_count > 0:
                    binned_positions[chrom].append(current_bin_start)
                    binned_coverages[chrom].append(bin_sum / bin_count)
                current_bin_start = pos
                current_bin_end = current_bin_start + bin_size
                bin_sum = cov
                bin_count = 1

        if bin_count > 0:
            binned_positions[chrom].append(current_bin_start)
            binned_coverages[chrom].append(bin_sum / bin_count)

    return binned_positions, binned_coverages

def generate_smoothed_coverage(bedgraph_file, bin_size=100000, sigma=50, mask_size=1000000):
    # Read the BEDGraph file
    positions, coverage_values = read_bedgraph(bedgraph_file)

    # Bin the data
    binned_positions, binned_coverages = bin_data(positions, coverage_values, bin_size, mask_size)

    # Concatenate positions and coverages for all chromosomes except Y
    concatenated_positions = []
    concatenated_coverages = []
    chromosome_starts = []
    current_position = 0
    for chrom in binned_positions.keys():
        if chrom == 'chrY':
            continue  # Skip chromosome Y
        
        chromosome_starts.append(current_position)
        concatenated_positions.extend([pos + current_position for pos in binned_positions[chrom]])
        concatenated_coverages.extend(binned_coverages[chrom])
        if binned_positions[chrom]:
            current_position += binned_positions[chrom][-1] + 1
    
    if concatenated_positions and concatenated_coverages:
        # Apply Gaussian smoothing with padding to reduce edge effects
        smoothed_coverage = smooth_with_padding(concatenated_coverages, sigma=sigma)
        return concatenated_positions, smoothed_coverage
    else:
        return None, None

# Function to apply Gaussian smoothing with padding to reduce edge effects
def smooth_with_padding(data, sigma):
    # Pad the data at both ends
    padded_data = np.pad(data, pad_width=sigma, mode='reflect')
    # Apply Gaussian smoothing
    smoothed_data = gaussian_filter1d(padded_data, sigma=sigma)
    # Remove the padding
    return smoothed_data[sigma:-sigma]


def get_max_y_ATAC_value(smoothed_coverage):
    if smoothed_coverage is None:
        return None
    
    max_y_value = np.max(smoothed_coverage)
    return max_y_value