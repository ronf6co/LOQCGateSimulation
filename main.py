from Cluster import Cluster
from FusionGate import FusionGate
from photonic_circuit import PhotonicCircuit
import numpy as np

if __name__ == '__main__':
    c1 = Cluster(2)
    c2 = Cluster(2)
    c3 = Cluster(2)

    fm = FusionGate(error_angle=-0.01, error_axis=[1, 0, 0])
    fp = FusionGate(error_angle=0.01, error_axis=[1, 0, 0])
    fi = FusionGate(error_angle=-0.0, error_axis=[1, 0, 0])

    co1pm = c1.fusion(other=c2, n1=1, n2=2, fusion_gate=f  m, clicks='00')
    co2pm = co1pm.fusion(other=c3, n1=1, n2=2, fusion_gate=fp, clicks='00')
    print(co2pm)

    co1pp = c1.fusion(other=c2, n1=1, n2=2, fusion_gate=fp, clicks='00')
    co2pp = co1pp.fusion(other=c3, n1=1, n2=2, fusion_gate=fp, clicks='00')
    print(co2pp)

    co1i = c1.fusion(other=c2, n1=1, n2=2, fusion_gate=fi, clicks='00')
    co2i = co1i.fusion(other=c3, n1=1, n2=2, fusion_gate=fi, clicks='00')
    print(co2i)

    print("PlusMinus dist from ideal : " + str(co2pm - co2i))
    print("PlusPlus dist from ideal : " + str(co2pp - co2i))