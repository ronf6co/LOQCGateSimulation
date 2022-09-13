import numpy as np
from dec2bin import *

from FusionGate import FusionGate


class Cluster:

    def __init__(self, n):
        self.res = 3
        if isinstance(n, int):
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

    def fusion(self, other, n1, n2, fusion_gate, clicks):
        # clicks - '00'/'01'/'10'/'11'
        # fusion_gate - object

        new_c_q_state = np.kron(self.q_state, other.q_state)

        ZZ = fusion_gate.U[int(clicks, 2)]

        U1 = np.eye(2 ** n1)
        U = np.kron(U1, ZZ)
        U2 = np.eye(2 ** (int(np.log2(len(new_c_q_state))) - n2 - 1))
        U = np.kron(U, U2)
        new_c_q_state = np.matmul(U, new_c_q_state)
        normalize = np.sqrt(sum(np.multiply(new_c_q_state, np.conj(new_c_q_state))))
        new_c_q_state = new_c_q_state / normalize

        phase1 = np.angle(new_c_q_state[0])
        new_c_q_state = new_c_q_state * np.exp(-1j * phase1)

        return Cluster(new_c_q_state)

    def __str__(self):
        s = ""
        s += "Size: " + str(int(np.log2(len(self.q_state))))
        s += "\nQ State:\n"
        for i in range(0, len(self.q_state)):
            s += str(round(self.q_state[i], self.res)) + "|" + format(i, '0' + str(
                int(np.log2(len(self.q_state)))) + 'b') + "⟩\n"
        return s

    def __sub__(self, other):

        # p1 = sum(np.multiply(self.q_state, np.conj(self.q_state)))
        p1 = np.outer(self.q_state, np.conj(self.q_state))
        # p2 = sum(np.multiply(other.q_state, np.conj(other.q_state)))
        p2 = np.outer(other.q_state, np.conj(other.q_state))

        return 0.5 * np.trace(abs(p1 - p2))



