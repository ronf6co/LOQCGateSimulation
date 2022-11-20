from Cluster import Cluster
from FusionGate import FusionGate

fi = FusionGate(error_angle=0, error_axis=[1,0,0])


for mini_c_size in [2]:
    for fusions_amount in [1]:

        c12i_00 = Cluster.fuseLinearClusters(size_of_cluster=mini_c_size,
                                          clicks="00",
                                          fusions=[fi] * fusions_amount,
                                          rot_state_angle=0,
                                          rot_state_axis=[1,0,0])

        c12i_01 = Cluster.fuseLinearClusters(size_of_cluster=mini_c_size,
                                             clicks="01",
                                             fusions=[fi] * fusions_amount,
                                             rot_state_angle=0,
                                             rot_state_axis=[1,0,0])

        c12i_ideal = Cluster(mini_c_size*(fusions_amount+1)-fusions_amount*2)

        print("Cluster of 2 is: \n" + str(c12i_ideal))
        print("Clister of 2 fusion with cluster of 2 with clicks 00: \n" + str(c12i_00))
        print("Clister of 2 fusion with cluster of 2 with clicks 01: \n" + str(c12i_01))