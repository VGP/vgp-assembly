from __future__ import print_function

import sys
from itertools import cycle

input_file = sys.argv[1]
number_of_split = int(sys.argv[2])

records = open(input_file).read().strip().split('\n')

output_name = input_file.replace('.bed', '').replace('.txt', '')
output_number = cycle(range(1, number_of_split+1))

for record in records:
    fdd = open(output_name + '_' + str(output_number.next()) + '.bed', 'a')
    fdd.writelines(record+'\n')
    fdd.close()
