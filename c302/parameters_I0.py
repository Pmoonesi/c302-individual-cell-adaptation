"""

Parameters C:
    Cells:           Single compartment, conductance based cell models with HH like ion channels
    Chem Synapses:   Event based, ohmic; one rise & one decay constant
    Gap junctions:   Electrical connection; current linerly depends on difference in voltages

ASSESSMENT:
    May be possible to use this to generate oscilliatory behaviour, but use of event based synapses normally requires
    spiking in cells, so core neurons will have to have clear spikes.

"""

from neuroml import GapJunction
from neuroml import PulseGenerator
from neuroml import FixedFactorConcentrationModel

from c302.bioparameters import c302ModelPrototype


class ParameterisedModel(c302ModelPrototype):
    def __init__(self):
        super(ParameterisedModel, self).__init__()
        self.level = "I"
        self.custom_component_types_definitions = ["interactome_cell.xml", "interactome_synapse.xml"]

        self.set_default_bioparameters()
        self.print_("Set default parameters for %s" % self.level)

    def set_default_bioparameters(self):

        self.add_bioparameter(
            "neuron_capacitance", "1.5 pF", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_leak_conductance", "10 pS", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_leak_potential", "-35 mV", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_activity_rise", "0.666 per_ms", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_activity_decay", "3.33 per_ms", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_sigmoid_beta", "0.125 per_mV", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_v_threshold", "0 mV", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "neuron_spike_threshold", "1000 mV", "BlindGuess", "0.1"
        )


        self.add_bioparameter(
            "neuron_to_neuron_chem_exc_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )
        self.add_bioparameter(
            "neuron_to_muscle_chem_exc_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )
        
        self.add_bioparameter("chem_exc_syn_erev", "0 mV", "BlindGuess", "0.1")
 
        self.add_bioparameter(
            "neuron_to_neuron_chem_inh_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )
        self.add_bioparameter(
            "neuron_to_muscle_chem_inh_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )

        self.add_bioparameter("chem_inh_syn_erev", "-48 mV", "BlindGuess", "0.1")

        self.add_bioparameter(
            "neuron_to_neuron_elec_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )
        self.add_bioparameter(
            "neuron_to_muscle_elec_syn_gbase", "100 pS", "BlindGuess", "0.1"
        )

        #NOTE - these are here only because the c302 requires a concentration model if your model is conductance based
        self.add_bioparameter(
            "ca_conc_decay_time", "13.811870945509265 ms", "BlindGuess", "0.1"
        )
        self.add_bioparameter(
            "ca_conc_rho", "0.000238919 mol_per_m_per_A_per_s", "BlindGuess", "0.1"
        )

        self.add_bioparameter(
            "ca_conc_decay_time_muscle", "60.811870945509265 ms", "BlindGuess", "0.1"
        )
        self.add_bioparameter(
            "ca_conc_rho_muscle",
            "0.002238919 mol_per_m_per_A_per_s",
            "BlindGuess",
            "0.1",
        )

        self.add_bioparameter(
            "unphysiological_offset_current", "0 pA", "KnownError", "0"
        )  # Can be activated later
        self.add_bioparameter(
            "unphysiological_offset_current_del", "0 ms", "KnownError", "0"
        )
        self.add_bioparameter(
            "unphysiological_offset_current_dur", "2000 ms", "KnownError", "0"
        )

    def create_models(self):
        self.create_generic_muscle_cell()
        self.create_generic_neuron_cell()
        self.create_offset_current()
        self.create_neuron_to_neuron_syn()
        self.create_neuron_to_muscle_syn()

    def create_generic_muscle_cell(self):
        
        C= self.get_bioparameter("neuron_capacitance").value
        G_c= self.get_bioparameter("neuron_leak_conductance").value
        E_cell= self.get_bioparameter("neuron_leak_potential").value
        a_r= self.get_bioparameter("neuron_activity_rise").value
        a_d= self.get_bioparameter("neuron_activity_decay").value
        beta= self.get_bioparameter("neuron_sigmoid_beta").value
        Vth= self.get_bioparameter("neuron_v_threshold").value
        spike_thresh= self.get_bioparameter("neuron_spike_threshold").value

        self.generic_muscle_cell = InteractomeNeuron(id="GenericMuscleCell",
                                                     C=C,
                                                     G_c=G_c,
                                                     E_cell=E_cell,
                                                     a_r=a_r,
                                                     a_d=a_d,
                                                     beta=beta,
                                                     Vth=Vth,
                                                     spike_thresh=spike_thresh)

    def create_generic_neuron_cell(self):
        C= self.get_bioparameter("neuron_capacitance").value
        G_c= self.get_bioparameter("neuron_leak_conductance").value
        E_cell= self.get_bioparameter("neuron_leak_potential").value
        a_r= self.get_bioparameter("neuron_activity_rise").value
        a_d= self.get_bioparameter("neuron_activity_decay").value
        beta= self.get_bioparameter("neuron_sigmoid_beta").value
        Vth= self.get_bioparameter("neuron_v_threshold").value
        spike_thresh= self.get_bioparameter("neuron_spike_threshold").value

        self.generic_neuron_cell = InteractomeNeuron(id="GenericNeuronCell",
                                                     C=C,
                                                     G_c=G_c,
                                                     E_cell=E_cell,
                                                     a_r=a_r,
                                                     a_d=a_d,
                                                     beta=beta,
                                                     Vth=Vth,
                                                     spike_thresh=spike_thresh)

    def create_offset_current(self):
        self.offset_current = PulseGenerator(
            id="offset_current",
            delay=self.get_bioparameter("unphysiological_offset_current_del").value,
            duration=self.get_bioparameter("unphysiological_offset_current_dur").value,
            amplitude=self.get_bioparameter("unphysiological_offset_current").value,
        )

    def create_neuron_to_neuron_syn(self):
        self.neuron_to_neuron_exc_syn = GradedSynapse3(
            id="neuron_to_neuron_exc_syn",
            conductance=self.get_bioparameter("neuron_to_neuron_chem_exc_syn_gbase").value,
            erev=self.get_bioparameter("chem_exc_syn_erev").value
        )

        self.neuron_to_neuron_inh_syn = GradedSynapse3(
            id="neuron_to_neuron_inh_syn",
            conductance=self.get_bioparameter("neuron_to_neuron_chem_inh_syn_gbase").value,
            erev=self.get_bioparameter("chem_inh_syn_erev").value
        )

        self.neuron_to_neuron_elec_syn = GapJunction(
            id="neuron_to_neuron_elec_syn",
            conductance=self.get_bioparameter("neuron_to_neuron_elec_syn_gbase").value,
        )

    def create_neuron_to_muscle_syn(self):
        self.neuron_to_muscle_exc_syn = GradedSynapse3(
            id="neuron_to_muscle_exc_syn",
            conductance=self.get_bioparameter("neuron_to_muscle_chem_exc_syn_gbase").value,
            erev=self.get_bioparameter("chem_exc_syn_erev").value
        )

        self.neuron_to_muscle_inh_syn = GradedSynapse3(
            id="neuron_to_muscle_inh_syn",
            conductance=self.get_bioparameter("neuron_to_muscle_chem_inh_syn_gbase").value,
            erev=self.get_bioparameter("chem_inh_syn_erev").value
        )

        self.neuron_to_muscle_elec_syn = GapJunction(
            id="neuron_to_muscle_elec_syn",
            conductance=self.get_bioparameter("neuron_to_muscle_elec_syn_gbase").value,
        )

    def get_cell(self, cell, type):
        self.found_specific_param = False

        C = self.get_cell_param(cell, "%s_%s", "neuron_%s", "capacitance")
        G_c = self.get_cell_param(cell, "%s_%s", "neuron_%s", "leak_conductance")
        E_cell = self.get_cell_param(cell, "%s_%s", "neuron_%s", "leak_potential")
        a_r = self.get_cell_param(cell, "%s_%s", "neuron_%s", "activity_rise")
        a_d = self.get_cell_param(cell, "%s_%s", "neuron_%s", "activity_decay")
        beta = self.get_cell_param(cell, "%s_%s", "neuron_%s", "sigmoid_beta")
        Vth = self.get_cell_param(cell, "%s_%s", "neuron_%s", "v_threshold")
        spike_thresh = self.get_cell_param(cell, "%s_%s", "neuron_%s", "spike_threshold")

        cell_obj = None

        if self.found_specific_param:
            cell_id = "%s_%s_cell" % (cell, type)
            cell_obj = InteractomeNeuron(id=cell_id,
                                C=C,
                                G_c=G_c,
                                E_cell=E_cell,
                                a_r=a_r,
                                a_d=a_d,
                                beta=beta,
                                Vth=Vth,
                                spike_thresh=spike_thresh)

        else:
            cell_obj = self.generic_neuron_cell if type == "neuron" else self.generic_muscle_cell
        
        return cell_obj

    def get_elec_syn(self, pre_cell, post_cell, type):
        self.found_specific_param = False
        if type == "neuron_to_neuron":
            gbase = self.get_conn_param(
                pre_cell, post_cell, "%s_to_%s_elec_syn_%s", "neuron_to_neuron_elec_syn_%s", "gbase",
            )
            conn_id = "neuron_to_neuron_elec_syn"
        elif type == "neuron_to_muscle":
            gbase = self.get_conn_param(
                pre_cell, post_cell, "%s_to_%s_elec_syn_%s", "neuron_to_muscle_elec_syn_%s", "gbase",
            )
            conn_id = "neuron_to_muscle_elec_syn"

        if self.found_specific_param:
            conn_id = "%s_to_%s_elec_syn" % (pre_cell, post_cell)

        return GapJunction(id=conn_id, conductance=gbase)

    def get_exc_syn(self, pre_cell, post_cell, type):
        self.found_specific_param = False

        specific_param_template = "%s_to_%s_chem_exc_syn_%s"
        if type == "neuron_to_neuron":
            gbase = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "neuron_to_neuron_chem_exc_syn_%s", "gbase",
            )
            erev = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "chem_exc_syn_%s", "erev"
            )

            conn_id = "neuron_to_neuron_exc_syn"

        elif type == "neuron_to_muscle":
            gbase = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "neuron_to_muscle_chem_exc_syn_%s", "gbase",
            )
            erev = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "chem_exc_syn_%s", "erev"
            )

            conn_id = "neuron_to_muscle_exc_syn"

        if self.found_specific_param:
            conn_id = "%s_to_%s_exc_syn" % (pre_cell, post_cell)

        return GradedSynapse3(
            id=conn_id, conductance=gbase, erev=erev
        )

    def get_inh_syn(self, pre_cell, post_cell, type):
        self.found_specific_param = False

        specific_param_template = "%s_to_%s_chem_inh_syn_%s"
        if type == "neuron_to_neuron":
            gbase = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "neuron_to_neuron_chem_inh_syn_%s", "gbase",
            )
            erev = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "chem_inh_syn_%s", "erev"
            )

            conn_id = "neuron_to_neuron_inh_syn"

        elif type == "neuron_to_muscle":
            gbase = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "neuron_to_muscle_chem_inh_syn_%s", "gbase",
            )
            erev = self.get_conn_param(
                pre_cell, post_cell, specific_param_template, "chem_inh_syn_%s", "erev"
            )

            conn_id = "neuron_to_muscle_inh_syn"

        if self.found_specific_param:
            conn_id = "%s_to_%s_inh_syn" % (pre_cell, post_cell)

        return GradedSynapse3(
            id=conn_id, conductance=gbase, erev=erev
        )

    def create_n_connection_synapse(self, prototype_syn, n, nml_doc, existing_synapses):
        if prototype_syn.id in existing_synapses:
            return existing_synapses[prototype_syn.id]

        existing_synapses[prototype_syn.id] = prototype_syn
        if isinstance(prototype_syn, GradedSynapse3):
            nml_doc.graded_synapses.append(prototype_syn)
        elif isinstance(prototype_syn, GapJunction):
            nml_doc.gap_junctions.append(prototype_syn)
        else:
            info = "%s conns; %s" % (n, existing_synapses)
            del existing_synapses[prototype_syn.id]
            raise Exception("Unknown synapse type: %s (%s)" % (prototype_syn.id, info))

        return prototype_syn

    def is_analog_conn(self, syn):
        return isinstance(syn, GradedSynapse3)

    def is_elec_conn(self, syn):
        return isinstance(syn, GapJunction)

class InteractomeNeuron:
    def __init__(self, id, C, G_c, E_cell, a_r, a_d, beta, Vth, spike_thresh):
        self.id = id
        self.C = C
        self.G_c = G_c
        self.E_cell = E_cell
        self.a_r = a_r
        self.a_d = a_d
        self.beta = beta
        self.Vth = Vth
        self.spike_thresh = spike_thresh

    def export(self, outfile, level, namespace, name_, pretty_print=True, **kwargs_):
        outfile.write(
            "    " * level
            + '<InteractomeNeuron id="%s" C="%s" G_c="%s" E_cell="%s" a_r="%s" a_d="%s" beta="%s" Vth="%s" spike_thresh="%s"/>\n'
            % (
                self.id,
                self.C,
                self.G_c,
                self.E_cell,
                self.a_r,
                self.a_d,
                self.beta,
                self.Vth,
                self.spike_thresh
            )
        )

class GradedSynapse3:
    def __init__(self, id, conductance, erev):
        self.id = id
        self.conductance = conductance
        self.erev = erev

    def export(self, outfile, level, namespace, name_, pretty_print=True, **kwargs_):
        outfile.write(
            "    " * level
            + '<GradedSynapse3 id="%s" conductance="%s" erev="%s"/>\n'
            % (
                self.id,
                self.conductance,
                self.erev,
            )
        )