SYMBOLIC = False
import sympy as sp
import numpy as np


class PhotonicCircuit:
    def __init__(self, n):
        self.n = n
        self.unitary = np.eye(self.n, dtype=np.complex128)

    def apply_coupler(self, coupler, wg1, wg2):
        if SYMBOLIC:
            coup_matrix = sp.eye(self.n, dtype=np.complex128)
        else:
            coup_matrix = np.eye(self.n, dtype=np.complex128)
        coup_matrix[wg1, wg1] = coupler[0, 0]
        coup_matrix[wg1, wg2] = coupler[0, 1]
        coup_matrix[wg2, wg1] = coupler[1, 0]
        coup_matrix[wg2, wg2] = coupler[1, 1]
        self.unitary = np.matmul(coup_matrix, self.unitary)

    def project_on_computation_basis(self, basis):
        m = len(basis)  # 2^Number of qubits

        if SYMBOLIC:
            logical_unitary = sp.zeros(m, m)
        else:
            logical_unitary = np.zeros((m, m), dtype=np.complex128)
        for i in range(m):
            for j in range(m):
                v_ind = [basis[i][0], basis[i][1]]
                u_ind = [basis[j][0], basis[j][1]]
                logical_unitary[i, j] = self.unitary[v_ind[0], u_ind[0]] * self.unitary[v_ind[1], u_ind[1]] \
                                        + self.unitary[v_ind[0], u_ind[1]] * self.unitary[v_ind[1], u_ind[0]]
        return logical_unitary
