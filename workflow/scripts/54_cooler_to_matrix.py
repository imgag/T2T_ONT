#!/usr/bin/env python3
# Converts a contact list (bin1 bin2 count) to a dense symmetric matrix for 3D structure tools

import sys
import numpy as np

# Read all contacts from stdin
contacts = []
max_bin = 0
min_bin = float('inf')

for line in sys.stdin:
    parts = line.strip().split()
    if len(parts) < 3:
        continue
    try:
        b1, b2, val = int(parts[0]), int(parts[1]), float(parts[2])
        contacts.append((b1, b2, val))
        max_bin = max(max_bin, b1, b2)
        min_bin = min(min_bin, b1, b2)
    except ValueError:
        continue

if not contacts:
    # Output a 1x1 zero matrix if no contacts
    np.savetxt(sys.stdout, np.zeros((1, 1)), fmt='%.6f', delimiter=' ')
    sys.exit(0)

# Determine matrix size based on max bin index (assuming 0-based bins)
n = max_bin - min_bin + 1
matrix = np.zeros((n, n))

for b1, b2, val in contacts:
    i = b1 - min_bin
    j = b2 - min_bin
    matrix[i, j] = val
    matrix[j, i] = val  # Symmetric

np.savetxt(sys.stdout, matrix, fmt='%.6f', delimiter=' ')