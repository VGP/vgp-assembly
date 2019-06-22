from __future__ import print_function

import csv
import sys

file_name = sys.argv[1]
fold_coverage = int(sys.argv[2])

with open(file_name, 'rd') as fd:
    fd_reader = csv.reader(fd)
    fd_columns = fd_reader.next()
    fd_value = fd_reader.next()
    fd_dict = dict(zip(fd_columns, fd_value))
    mean_depth = fd_dict['mean_depth']

print(int(round(float(mean_depth)*fold_coverage)))
