from photonic_circuit import PhotonicCircuit
import numpy as np

if __name__ == '__main__':
    I = np.eye(2)
    X = np.array([[0, 1], [1, 0]])
    X_sqrt = np.cos(np.pi/4) * I - np.sin(np.pi/4) * 1j * X
    circuit = PhotonicCircuit(4)
    circuit.apply_coupler(X_sqrt, 0, 1)
    circuit.apply_coupler(X_sqrt, 2, 3)
    circuit.apply_coupler(-1j*X, 1, 2)
    circuit.apply_coupler(X_sqrt, 0, 1)
    circuit.apply_coupler(X_sqrt, 2, 3)
    basis = [[0, 2], [0, 3], [1, 2], [1, 3]]
    U = circuit.project_on_computation_basis(basis)
    print(U)