from __future__ import print_function

import sys

filename = sys.argv[1]
extension = sys.argv[2]
suffix = '.'.join(filename.split('.')[1:])
print(filename.split('_')[0]+extension+'.'+suffix)