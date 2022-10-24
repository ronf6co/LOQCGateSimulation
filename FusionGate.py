SYMBOLIC = True

import sympy as sp
import numpy as np

from photonic_circuit import PhotonicCircuit

if SYMBOLIC:
    I = sp.eye(2)
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -1j], [1j, 0]])
    Z = sp.Matrix([[1, 0], [0, -1]])
else:
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.array([[1, 0], [0, -1]])


class FusionGate:

    def __init__(self, error_angle=0, error_axis=None):
        if SYMBOLIC:
            if error_axis is None:
                error_axis = [1, 0, 0]

            self.error_angle = error_angle

            self.R_error = sp.cos(error_angle / 2) * I + 1j * sp.sin(error_angle / 2) * \
                      (error_axis[0] * X + error_axis[0] * Y + error_axis[0] * Z)
            X_sqrt = sp.cos(np.pi / 4) * I - sp.sin(np.pi / 4) * 1j * X

            circuit = PhotonicCircuit(4)
            circuit.apply_coupler(X_sqrt * self.R_error, 0, 1)
            circuit.apply_coupler(X_sqrt * self.R_error, 2, 3)
            circuit.apply_coupler(-1j * X_sqrt * self.R_error, 1, 2)
            circuit.apply_coupler(X_sqrt * self.R_error, 0, 1)
            circuit.apply_coupler(X_sqrt * self.R_error, 2, 3)
            self.basis = [[0, 2], [0, 3], [1, 2], [1, 3]]
            self.U = circuit.project_on_computation_basis(self.basis)

        else:
            if error_axis is None:
                error_axis = [1, 0, 0]

            self.error_angle = error_angle

            self.R_error = np.cos(error_angle / 2) * I + 1j * np.sin(error_angle / 2) * \
                           (error_axis[0] * X + error_axis[0] * Y + error_axis[0] * Z)
            X_sqrt = np.cos(np.pi / 4) * I - np.sin(np.pi / 4) * 1j * X

            circuit = PhotonicCircuit(4)
            circuit.apply_coupler(np.matmul(X_sqrt, self.R_error), 0, 1)
            circuit.apply_coupler(np.matmul(X_sqrt, self.R_error), 2, 3)
            circuit.apply_coupler(-1j * np.matmul(X, self.R_error), 1, 2)
            circuit.apply_coupler(np.matmul(X_sqrt, self.R_error), 0, 1)
            circuit.apply_coupler(np.matmul(X_sqrt, self.R_error), 2, 3)
            self.basis = [[0, 2], [0, 3], [1, 2], [1, 3]]
            self.U = circuit.project_on_computation_basis(self.basis)

    def __str__(self):
        s = ""
        s += "\n______________________Fusion:__________________\n"
        s += "U : \n"
        if self.error_angle == 0 or SYMBOLIC:
            round_res = 7
        else:
            round_res = int(7+np.log10(1/self.error_angle))
        if not SYMBOLIC:
            s += str(np.round(self.U, round_res))
        else:
            s += str(sp.simplify(self.U))
        s += "\n\nBasis :\n"
        s += str(self.basis)

        s += "\n\nR error : \n"
        if not SYMBOLIC:
            s += str(np.round(self.R_error, round_res))
        else:
            s += str(self.R_error)
        s += "\n_________________________________________________"
        return s


# print(FusionGate(error_angle=0.01, error_axis=[0, 0, 1]))
