
import scvelo as scv
import sys

loom = sys.argv[1]
test = scv.read(loom)
test.write('%s.h5ad' % loom)

