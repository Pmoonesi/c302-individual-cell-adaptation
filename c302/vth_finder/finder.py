import numpy as np
from scipy import sparse, linalg

import os
import json

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
def get_path(relative):
    return os.path.join(BASE_DIR, relative)

class Integrator:

    def __init__(self):

        """ Number of Neurons """
        self.N = 279

        """ Cell membrane conductance (pS) """
        self.Gc = 0.1 # again this is divided by 100pS so actually this is 10pS

        """ Cell Membrane Capacitance """ # should be 1 pF so this is again 100 times less than expected. (if only the current is also 100 times less than expected, they cancel each other out and all is good)
        # but this is divided by 100pS so instead of 1.5pF we have 0.015 arb
        self.C = 0.015

        """ Gap Junctions (Electrical, 279*279) """ # These stand for the conductances. They are counts.
        # this is divided by 100pS so this is 100pS in reality
        self.ggap = 1.0 # The coefficients aren't exact but the ratio's make sense. (gap junctions and synapses should be around 100 pS and the cell membranes 10 pS).
        self.Gg_Static = np.load(get_path('Gg.npy'))

        """ Synaptic connections (Chemical, 279*279) """ # These stand for the conductances. They are counts.
        self.gsyn = 1.0
        self.Gs_Static = np.load(get_path('Gs.npy'))

        """ Leakage potential (mV) """
        self.Ec = -35.0

        """ Directionality (279*1) """ # the threshold/potential in synapses (why only 26 non zero?)
        self.E = np.load(get_path('emask.npy'))
        self.E = -48.0 * self.E
        self.EMat = np.tile(np.reshape(self.E, self.N), (self.N, 1))

        """ Synaptic Activity Parameters """
        self.ar = 1.0/1.5 # Synaptic activity's rise time
        self.ad = 5.0/1.5 # Synaptic activity's decay time
        self.B = 0.125 # Width of the sigmoid (mv^-1)

        """ Input_Mask/Continuous Transtion """
        self.transit_Mat = np.zeros((2, self.N))

        self.t_Tracker = 0 # don't know what this is, suppose this keeps track of time
        self.Iext = 100000  # Our input is 10^-15 A, we multiply it by Iext to get 10^-10 and then we have our number that is divided by 100 like the rest.

        self.rate = 0.025 # the denominator
        self.offset = 0.15 # the offset of transition after switch is 

        """ Connectome Arrays """
        self.Gg_Dynamic = self.Gg_Static.copy()
        self.Gs_Dynamic = self.Gs_Static.copy()

        self.EffVth(self.Gg_Static, self.Gs_Static)

        self.load_presets()
    
    def load_presets(self):
        # Load all presets needed
        self.presets = {'RST': self.load_preset('RST')}

    def load_preset(self, preset_name):
        full_path = get_path(os.path.join('presets', preset_name + '.json'))
        if not os.path.isfile(full_path):
            return {}
        try:
            with open(full_path, 'rb') as file:
                return json.load(file)
        except Exception as e:
            print(f"Error loading preset {preset_name}: {e}")
            return {}

    """ Efficient V-threshold computation """
    def EffVth(self, Gg, Gs):

        Gcmat = np.multiply(self.Gc, np.eye(self.N))
        EcVec = np.multiply(self.Ec, np.ones((self.N, 1)))

        M1 = -Gcmat
        b1 = np.multiply(self.Gc, EcVec)

        Ggap = np.multiply(self.ggap, Gg)
        Ggapdiag = np.subtract(Ggap, np.diag(np.diag(Ggap)))
        Ggapsum = Ggapdiag.sum(axis = 1)
        Ggapsummat = sparse.spdiags(Ggapsum, 0, self.N, self.N).toarray()
        M2 = -np.subtract(Ggapsummat, Ggapdiag)

        Gs_ij = np.multiply(self.gsyn, Gs)
        s_eq = round((self.ar/(self.ar + 2 * self.ad)), 4)
        sjmat = np.multiply(s_eq, np.ones((self.N, self.N)))
        S_eq = np.multiply(s_eq, np.ones((self.N, 1))) #?
        Gsyn = np.multiply(sjmat, Gs_ij)
        Gsyndiag = np.subtract(Gsyn, np.diag(np.diag(Gsyn)))
        Gsynsum = Gsyndiag.sum(axis = 1)
        M3 = -sparse.spdiags(Gsynsum, 0, self.N, self.N).toarray()

        b3 = np.dot(Gs_ij, np.multiply(s_eq, self.E))

        M = M1 + M2 + M3

        (self.P, self.LL, self.UU) = linalg.lu(M) # this part is always the same so it does some precomputing (i guess lu decomposition)

        bbb = -b1 - b3 # once we add the external current to this we can solve for Vth
        self.bb = np.reshape(bbb, self.N)

    def EffVth_rhs(self, Iext, InMask):

        InputMask = np.multiply(Iext, InMask)
        self.b = np.subtract(self.bb, InputMask) # this is where we complete the b in our Ax = b and then we solve for x (Vth).

        Vth = linalg.solve_triangular(self.UU, linalg.solve_triangular(self.LL, self.b, lower = True, check_finite=False), check_finite=False)

        return Vth

    def solve(self, stimuli, names=False):

        vth = self.EffVth_rhs(self.Iext, stimuli)

        if names == False:

            return vth
        
        neuron_names = np.loadtxt(get_path("neuron_names.npy"), dtype=str)

        return dict(zip(neuron_names, vth))

    def evaluate(self, x):

        A=self.P @ self.LL @ self.UU

        temp = A @ x - self.b

        return np.linalg.norm(temp)

    def get_stimuli(self, stim_object):

        self.load_presets()
        preset = self.presets['RST']
        for neuron, current in stim_object.items():
            preset[neuron]['inputCurrent'] = current
        stimuli_currents = [entry["inputCurrent"] for entry in preset.values()]

        return stimuli_currents

    def get_seq(self):

        s_eq = round((self.ar/(self.ar + 2 * self.ad)), 4)
        return np.multiply(s_eq, np.ones((self.N, 1)))

    def compute_jacobian(self, y, stimuli):

        Vvec, SVec = np.split(y, 2)
        Vrep = np.tile(Vvec, (self.N, 1))

        J1_M1 = -np.multiply(self.Gc, np.eye(self.N))
        Ggap = np.multiply(self.ggap, self.Gg_Dynamic)
        Ggapsumdiag = -np.diag(Ggap.sum(axis = 1))
        J1_M2 = np.add(Ggap, Ggapsumdiag) 
        Gsyn = np.multiply(self.gsyn, self.Gs_Dynamic)
        J1_M3 = np.diag(np.dot(-Gsyn, SVec))

        J1 = (J1_M1 + J1_M2 + J1_M3) / self.C

        J2_M4_2 = np.subtract(self.EMat, np.transpose(Vrep))
        J2 = np.multiply(Gsyn, J2_M4_2) / self.C

        Vth = self.EffVth_rhs(self.Iext, stimuli)

        ## overflow: exp goes to infinite and the expression goes to zero anyways.
        try:
            sigmoid_V = np.reciprocal(1.0 + np.exp(-self.B*(np.subtract(Vvec, Vth))))
        except FloatingPointError as fpe:
            # this is just a reminder that the overflow problem exists
            sigmoid_V = np.zeros(279)
        
        J3_1 = np.multiply(self.ar, 1 - SVec)
        J3_2 = np.multiply(self.B, sigmoid_V)
        J3_3 = 1 - sigmoid_V
        J3 = np.diag(np.multiply(np.multiply(J3_1, J3_2), J3_3))

        J4 = np.diag(np.subtract(np.multiply(-self.ar, sigmoid_V), self.ad))

        J_row1 = np.hstack((J1, J2))
        J_row2 = np.hstack((J3, J4))
        J = np.vstack((J_row1, J_row2))

        return J

    def is_steq_stable(self, stimuli):
        vth = self.solve(stimuli)
        seq = self.get_seq()
        y = np.vstack([vth.reshape(-1,1), seq]).flatten()
        J = self.compute_jacobian(y, stimuli)
        eigenvalues, eigenvectors = np.linalg.eig(J)
        real_eigs = np.array(list(map(lambda x: x.real, eigenvalues)))
        all_real_parts_negative = np.all(real_eigs <= 0)
        return all_real_parts_negative

if __name__ == "__main__":

    integrator = Integrator()
    stim_obj = {"PLMR": 5e-4, "PLML": 5e-4}
    stim = integrator.get_stimuli(stim_obj)
    vth = integrator.solve(stim)

    print(vth)

