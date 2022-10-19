import numpy as np

import Cluster


class Star4Cluster:

    def __init__(self, n):
        #   4
        # 1 5 3
        #   2
        self.res = 3
        x = (1 / 2 ** 0.5) * np.array([1, 1])
        for _ in range(1, 5):
            x = np.kron(x, ([1, 1]))

        x = Cluster.CZ(x, 1, 5)
        x = Cluster.CZ(x, 2, 5)
        x = Cluster.CZ(x, 3, 5)
        x = Cluster.CZ(x, 4, 5)

        # normalize = np.sqrt(sum(x * np.conj(x)))
        # x = x / normalize

        self.q_state = x

    def fusion(self, other, fusion_gate, clicks):
        # clicks - '00'/'01'/'10'/'11'
        # fusion_gate - object

        zero @



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

