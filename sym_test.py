SYMBOLIC = True

import sympy as sp
import numpy as np

from Cluster import Cluster
from FusionGate import FusionGate
from photonic_circuit import PhotonicCircuit


error_angle = sp.symbols('t')
error_axis = [1, 0, 0]
fm = FusionGate(error_angle=-error_angle, error_axis=error_axis)
fp = FusionGate(error_angle=error_angle, error_axis=error_axis)
fi = FusionGate(error_angle=0, error_axis=error_axis)
print("Fusions:")
print(fm, fi, fp)

print("Clusters:")
# c12p = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fp])
# print(c12p)
# c12m = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fm])
# print(c12m)
# c12i = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fi])
# print(c12i)




# Worst +-+ Best -++
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='01', fusions_amount=3)
#
c12naive = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fp,fm,fp])
print(c12naive)
c12pm = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fm, fp, fp])
print(c12pm)
c12i = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fi]*3)
print(c12i)
print("PM dist from ideal : " + str(c12pm - c12i))
print("Naive dist from ideal : " + str(c12naive - c12i))