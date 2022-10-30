SYMBOLIC = False
import sympy as sp
import numpy as np

from Cluster import Cluster
from FusionGate import FusionGate
from photonic_circuit import PhotonicCircuit

def best_combo_with_vs_without_rotation(fi,fp,fm,cluster_size,clicks,fusions_amount, rot_state_angle, rot_state_axis=[1,0,0]):
    runOverPMCombinations(fi, fp, fm, cluster_size, clicks, fusions_amount, rot_state_angle=0, rot_state_axis=rot_state_axis)
    runOverPMCombinations(fi, fp, fm, cluster_size, clicks, fusions_amount, rot_state_angle=rot_state_angle, rot_state_axis=rot_state_axis)

def runOverPMCombinations(fi,fp,fm,cluster_size,clicks,fusions_amount, rot_state_angle=0, rot_state_axis=[1,0,0]):
    print("---------------------------------------")
    print("Finding combination for ",str(fusions_amount+1),"cluster of",str(cluster_size),"  Clicks are ",clicks)
    print("State change from cluster is theta=", str(rot_state_angle))
    print("Statrting...")
    c12naive = Cluster.fuseLinearClusters(size_of_cluster=cluster_size, clicks=clicks, fusions=[fp]*fusions_amount, rot_state_angle=rot_state_angle, rot_state_axis=rot_state_axis)
    c12i = Cluster.fuseLinearClusters(size_of_cluster=cluster_size, clicks=clicks, fusions=[fi] * fusions_amount, rot_state_angle=rot_state_angle, rot_state_axis=rot_state_axis)
    worst_case = float("-Inf")
    best_case = float("Inf")
    for i in range(2**fusions_amount):
        b = format(i, '0'+str(fusions_amount)+'b')
        fusions_list = []
        for bit in b:
            if bit == '0':
                fusions_list += [fp]
            else:
                fusions_list += [fm]
        c12pm = Cluster.fuseLinearClusters(size_of_cluster=cluster_size, clicks=clicks, fusions=fusions_list, rot_state_angle=rot_state_angle, rot_state_axis=rot_state_axis)
        # print("for: "+str(b)+" PM dist from ideal : " + str(c12pm - c12i))
        if c12pm - c12i > worst_case:
            worst_case = c12pm - c12i
            worst_case_raw = b
        if c12pm - c12i < best_case:
            best_case = c12pm - c12i
            best_case_raw = b

    print("Naive (00..0) dist from ideal : " + str(c12naive - c12i))

    print("Worst case is : " + worst_case_raw + " value: " + str(worst_case))
    print("Best case is : " + best_case_raw + " value: " + str(best_case))
    print("---------------------------------------")

if __name__ == '__main__':
    print("Init Fusion Gates...")
# a
    error_angle = 0.1
    error_axis = [0, 1, 0]
    fm = FusionGate(error_angle=-error_angle, error_axis=error_axis)
    fp = FusionGate(error_angle=error_angle, error_axis=error_axis)
    fi = FusionGate(error_angle=0, error_axis=error_axis)

    # print("O-O ~ O-O")
    # c12p = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fp])
    # print("O-O ~ O-O")
    # c12m = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fm])
    # c12i = Cluster.fuseLinearClusters(size_of_cluster=2, clicks='00', fusions=[fi])
    # print("O-O ~ O-O")
    # c12naive = Cluster.fuseLinearClusters(size_of_cluster=3, clicks='00', fusions=[fp,fp,fp,fp,fp])
    # c12pm = Cluster.fuseLinearClusters(size_of_cluster=3, clicks='00', fusions=[fp, fm, fp, fm, fp])
    # c12i = Cluster.fuseLinearClusters(size_of_cluster=3, clicks='00', fusions=[fi]*5)
    # print("PM dist from ideal : " + str(c12pm - c12i))
    # print("Naive dist from ideal : " + str(c12naive - c12i))

    # Like
    # Batch - Cluster(2) Clicks 00/11
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=3, rot_state_angle=0, rot_state_axis=[1,0,0])
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=6, rot_state_angle=0, rot_state_axis=[1,0,0])
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=8, rot_state_angle=0, rot_state_axis=[1,0,0])
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='11', fusions_amount=8, rot_state_angle=0, rot_state_axis=[1,0,0])

    best_combo_with_vs_without_rotation(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=3, rot_state_angle=0.2,
                          rot_state_axis=[1, 0, 0])
    best_combo_with_vs_without_rotation(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=6, rot_state_angle=0.2,
                          rot_state_axis=[1, 0, 0])
    best_combo_with_vs_without_rotation(fi, fp, fm, cluster_size=2, clicks='00', fusions_amount=8, rot_state_angle=0.2,
                          rot_state_axis=[1, 0, 0])
    best_combo_with_vs_without_rotation(fi, fp, fm, cluster_size=2, clicks='11', fusions_amount=8, rot_state_angle=0.2,
                          rot_state_axis=[1, 0, 0])

    # Batch - Cluster(2) Clicks 10
    # Worst ++++++ Best +-++---+
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='10', fusions_amount=8)
    # Worst ----- Best --+++
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='10', fusions_amount=5)
    # Worst --- Best ++-
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='10', fusions_amount=3)

    # Batch - Cluster(2) Clicks 01 odd
    # Worst +-+ Best -++
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='01', fusions_amount=3)
    # Worst +-+-+ Best ---++
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='01', fusions_amount=5)

    # Batch - Cluster(2) Clicks 01 Even
    # Worst ++++ Best -+-+
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='01', fusions_amount=4)
    # Worst ++++++ Best ++--+-
    # runOverPMCombinations(fi, fp, fm, cluster_size=2, clicks='01', fusions_amount=6)

    # Batch - Cluster(3) Clicks 00 Even/Odd
    # Worst ------, Best ++++++
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='00', fusions_amount=6)
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='00', fusions_amount=4)
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='00', fusions_amount=5)
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='00', fusions_amount=7)

    # Batch - Cluster(3) Clicks 11 Even/Odd
    # Worst ---, Best +++
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='11', fusions_amount=4)

    # Batch - Cluster(3) Clicks 01 Even/Odd (no big gap)
    # Worst +++++, Best -----
    # runOverPMCombinations(fi, fp, fm, cluster_size=3, clicks='01', fusions_amount=5)






