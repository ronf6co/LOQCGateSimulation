SYMBOLIC = False
import numpy as np

from Cluster import Cluster
from FusionGate import FusionGate
from photonic_circuit import PhotonicCircuit

from main import runOverPMCombinations

print("Init Fusion Gates...")
error_angle = 0.05
error_axis = [1, 0, 0]
fm = FusionGate(error_angle=-error_angle, error_axis=error_axis)
fp = FusionGate(error_angle=error_angle, error_axis=error_axis)
fi = FusionGate(error_angle=0, error_axis=error_axis)

