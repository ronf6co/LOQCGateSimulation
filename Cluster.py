SYMBOLIC = False
import sympy as sp
from sympy.physics.quantum import TensorProduct
import numpy as np
from dec2bin import *
from sympy.core.rules import Transform

from FusionGate import FusionGate

I = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])


class Cluster:

    def __init__(self, n, U=None):
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
        if U is not None:
            self.U = U
        else:
            self.U = np.eye(len(self.q_state))

    @staticmethod
    def CZ(x, n1, n2):
        for i in range(len(x)):
            b = format(i, '0' + str(int(np.log2(len(x)))) + 'b')
            if b[n1] == '1' and b[n2] == '1':
                x[i] = -x[i]
        return x

    def rotate(self, index, theta, axis):
        matrix1 = np.eye(2 ** index)
        matrix2 = np.cos(theta / 2) * I + 1j * np.sin(theta / 2) * \
                  (axis[0] * X + axis[1] * Y + axis[2] * Z)
        matrix3 = np.eye(2 ** (int(np.log2(len(self.q_state))) - index - 1))
        matrix_total = np.kron(matrix1, np.kron(matrix2, matrix3))
        self.q_state = np.matmul(matrix_total, self.q_state)

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
            normalize = np.sqrt(
                sum(np.multiply(new_c_q_state, np.conj(new_c_q_state))))
            U = U * (1 / normalize)
            new_c_q_state = U * new_c_q_state

            # normalize = sp.sqrt(sum(sp.matrices.dense.matrix_multiply_elementwise(new_c_q_state, sp.conjugate(new_c_q_state))))
            # new_c_q_state = new_c_q_state / normalize

            # phase1 = sp.log(new_c_q_state[0]).as_real_imag()[1]
            # new_c_q_state = new_c_q_state * sp.exp(-1j * phase1)

            # TODO U for symbolic
            raise NotImplementedError("SYMBOLIC not implemented")
        else:

            n1 = int(np.log2(len(self.q_state))) - 1
            U1 = np.eye(2 ** n1)
            n2 = int(np.log2(len(other.q_state))) - 1
            U2 = np.eye(2 ** n2)

            new_c_q_state = np.kron(self.q_state, other.q_state)

            ZZ = fusion_gate.U[int(clicks, 2)]

            U = np.kron(U1, ZZ)
            U = np.kron(U, U2)
            normalize = np.sqrt(
                sum(np.multiply(new_c_q_state, np.conj(new_c_q_state))))
            U = U * (1 / normalize)
            new_c_q_state = np.matmul(U, new_c_q_state)

            phase1 = np.angle(new_c_q_state[0])
            new_c_q_state = new_c_q_state * np.exp(-1j * phase1)

            original_u = np.kron(self.U, other.U)
            mid_u = np.kron(np.kron(np.eye(int(len(self.q_state) * 0.5)), ZZ),
                            np.eye(int(len(other.q_state) * 0.5)))
            new_U = np.matmul(mid_u, original_u)

        return Cluster(new_c_q_state, new_U)

    @staticmethod
    def fuseLinearClusters(size_of_cluster, clicks, fusions, rot_state_angle=0,
                           rot_state_axis=[1, 0, 0]):
        if isinstance(size_of_cluster, list):
            # assert(len(size_of_cluster) == len(fusions)+1)
            c_out = Cluster(size_of_cluster[0])
            c_out.rotate(1, rot_state_angle, rot_state_axis)
            for i in range(1, len(fusions) + 1):
                other_cluster = Cluster(size_of_cluster[i])
                other_cluster.rotate(1, rot_state_angle, rot_state_axis)
                c_out = c_out.fusion(other=other_cluster,
                                     fusion_gate=fusions[i - 1], clicks=clicks)
            return c_out
        else:
            c_out = Cluster(size_of_cluster)
            c_out.rotate(1, rot_state_angle, rot_state_axis)
            for i in range(1, len(fusions) + 1):
                other_cluster = Cluster(size_of_cluster)
                other_cluster.rotate(1, rot_state_angle, rot_state_axis)
                c_out = c_out.fusion(other=Cluster(size_of_cluster),
                                     fusion_gate=fusions[i - 1], clicks=clicks)
            return c_out

    def __str__(self):
        s = ""
        s += "Size: " + str(int(np.log2(len(self.q_state))))
        s += "\nQ State:\n"
        for i in range(0, len(self.q_state)):
            if SYMBOLIC:
                x = sp.simplify(self.q_state[i].xreplace(
                    Transform(lambda x: x.round(self.res),
                              lambda x: isinstance(x, float))))
                s += str(x) + "|" + format(i, '0' + str(
                    int(np.log2(len(self.q_state)))) + 'b') + "⟩\n"
            else:
                s += str(round(self.q_state[i], self.res)) \
                     + "|" + format(i,
                                    '0' + str(
                                        int(np.log2(
                                            len(self.q_state)))) + 'b') + "⟩\n"
        return s

    def __sub__(self, other):

        diff_mat = self.U - other.U
        return np.abs(np.sqrt(
            np.trace(np.matmul(diff_mat, np.conj(np.transpose(diff_mat))))))

        # yaron = (1/len(self.U))*np.trace(np.matmul(np.conj(np.transpose(self.U)), (self.U-other.U)))
        # return yaron

        # if SYMBOLIC:
        #     p1 = sp.outer(self.q_state, sp.conjugate(self.q_state))
        #     p2 = sp.outer(other.q_state, sp.conjugate(other.q_state))
        #
        #     return sp.simplify(0.5 * sp.trace(abs(p1 - p2)))
        # else:
        #
        #     # p1 = sum(np.multiply(self.q_state, np.conj(self.q_state)))
        #     p1 = np.outer(self.q_state, np.conj(self.q_state))
        #     # p2 = sum(np.multiply(other.q_state, np.conj(other.q_state)))
        #     p2 = np.outer(other.q_state, np.conj(other.q_state))
        #
        #     return 0.5 * np.trace(abs(p1 - p2))
