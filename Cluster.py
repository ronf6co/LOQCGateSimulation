SYMBOLIC = True
import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
from dec2bin import *
from sympy.core.rules import Transform

from FusionGate import FusionGate


class Cluster:

    def __init__(self, n):
        self.res = 3
        if isinstance(n, int):
            if SYMBOLIC:
                x = (1 / 2 ** 0.5) * sp.Matrix([1, 1])
                for _ in range(1, n):
                    x = TensorProduct(x, sp.Matrix([1, 1]))

                for i in range(0, n - 1):
                    x = Cluster.CZ(x, i, i + 1)

                # normalize = np.sqrt(sum(x * np.conj(x)))
                # x = x / normalize

                self.q_state = x
            else:
                x = (1 / 2 ** 0.5) * np.array([1, 1])
                for _ in range(1, n):
                    x = np.kron(x, ([1, 1]))

                for i in range(0, n - 1):
                    x = Cluster.CZ(x, i, i + 1)

                # normalize = np.sqrt(sum(x * np.conj(x)))
                # x = x / normalize

                self.q_state = x
        else:
            self.q_state = n

    @staticmethod
    def CZ(x, n1, n2):
        for i in range(len(x)):
            b = format(i, '0' + str(int(np.log2(len(x)))) + 'b')
            if b[n1] == '1' and b[n2] == '1':
                x[i] = -x[i]
        return x

    def fusion(self, other, fusion_gate, clicks):
        # clicks - '00'/'01'/'10'/'11'
        # fusion_gate - object

        if SYMBOLIC:
            n1 = int(np.log2(len(self.q_state))) - 1
            U1 = sp.eye(2 ** n1)
            n2 = int(np.log2(len(other.q_state))) - 1
            U2 = sp.eye(2 ** n2)

            new_c_q_state = TensorProduct(self.q_state, other.q_state)

            # ZZ = fusion_gate.U[int(clicks, 2)]
            ZZ = fusion_gate.U.row(int(clicks, 2))

            U = TensorProduct(U1, sp.Matrix([ZZ]))
            U = TensorProduct(U, U2)
            new_c_q_state = U * new_c_q_state


            # normalize = sp.sqrt(sum(sp.matrices.dense.matrix_multiply_elementwise(new_c_q_state, sp.conjugate(new_c_q_state))))
            # new_c_q_state = new_c_q_state / normalize

            # phase1 = sp.log(new_c_q_state[0]).as_real_imag()[1]
            # new_c_q_state = new_c_q_state * sp.exp(-1j * phase1)
        else:

            n1 = int(np.log2(len(self.q_state)))-1
            U1 = np.eye(2 ** n1)
            n2 = int(np.log2(len(other.q_state)))-1
            U2 = np.eye(2 ** n2)

            new_c_q_state = np.kron(self.q_state, other.q_state)

            ZZ = fusion_gate.U[int(clicks, 2)]

            U = np.kron(U1, ZZ)
            U = np.kron(U, U2)
            new_c_q_state = np.matmul(U, new_c_q_state)
            normalize = np.sqrt(sum(np.multiply(new_c_q_state, np.conj(new_c_q_state))))
            new_c_q_state = new_c_q_state / normalize

            phase1 = np.angle(new_c_q_state[0])
            new_c_q_state = new_c_q_state * np.exp(-1j * phase1)

        return Cluster(new_c_q_state)

    @staticmethod
    def fuseLinearClusters(size_of_cluster, clicks, fusions):
        c_out = Cluster(size_of_cluster)
        for i in range(1, len(fusions)+1):
            c_out = c_out.fusion(other=Cluster(size_of_cluster), fusion_gate=fusions[i - 1], clicks=clicks)
        return c_out

    def __str__(self):
        s = ""
        s += "Size: " + str(int(np.log2(len(self.q_state))))
        s += "\nQ State:\n"
        for i in range(0, len(self.q_state)):
            if SYMBOLIC:
                x = sp.simplify(self.q_state[i].xreplace(Transform(lambda x: x.round(self.res), lambda x: isinstance(x, float))))
                s += str(x) + "|" + format(i, '0' + str(
                    int(np.log2(len(self.q_state)))) + 'b') + "⟩\n"
            else:
                s += str(round(self.q_state[i], self.res)) + "|" + format(i, '0' + str(
                    int(np.log2(len(self.q_state)))) + 'b') + "⟩\n"
        return s

    def __sub__(self, other):

        if SYMBOLIC:
            p1 = sp.outer(self.q_state, sp.conjugate(self.q_state))
            p2 = sp.outer(other.q_state, sp.conjugate(other.q_state))

            return sp.simplify(0.5 * sp.trace(abs(p1 - p2)))
        else:

            # p1 = sum(np.multiply(self.q_state, np.conj(self.q_state)))
            p1 = np.outer(self.q_state, np.conj(self.q_state))
            # p2 = sum(np.multiply(other.q_state, np.conj(other.q_state)))
            p2 = np.outer(other.q_state, np.conj(other.q_state))

            return 0.5 * np.trace(abs(p1 - p2))



