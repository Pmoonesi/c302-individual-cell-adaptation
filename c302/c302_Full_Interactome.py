import c302
import sys
from c302.vth_finder.finder import Integrator
from c302.bioparameters import split_neuroml_quantity
from c302 import get_cell_names_and_connection

import random

import importlib
import re

neuron_pattern = re.compile(r'([A-Z]+)(\d*)')

def get_c302_name(neuron_name):
    match_result = neuron_pattern.search(neuron_name)
    name, number = match_result.group(1), match_result.group(2)

    if not number:
        return neuron_name
    
    number = int(number)

    return "%s%s" % (name, number)


def set_vth_bioparams(params, cells_to_stimulate, offset_current):

    magnitude, unit = split_neuroml_quantity(offset_current)

    stimulus_size = magnitude 

    if unit == "nA":
        stimulus_size *= 0.1
    elif unit == "pA":
        stimulus_size *= 0.0001

    stimulus_array = [stimulus_size] * len(cells_to_stimulate)

    stim_obj = dict(zip(cells_to_stimulate, stimulus_array))

    integrator = Integrator()
    stim = integrator.get_stimuli(stim_obj)
    vth_dict = integrator.solve(stim, names=True)

    for neuron, vth in vth_dict.items():
        params.add_bioparameter(
            "%s_v_threshold" % (get_c302_name(neuron)), "%.4f mV" % (vth),  "Testing Full net", "0"
        )

def set_init_potential_bioparams(params, base, data_reader="KimDataReaderV3"):
    magnitude, unit = split_neuroml_quantity(base)
    names, conns = get_cell_names_and_connection(data_reader)

    for neuron in names:

        randomized_leak_pot = magnitude + random.normalvariate(0, 1) * 1e-4

        params.add_bioparameter(
            "%s_v_init" % (get_c302_name(neuron)), "%.6f mV" % (randomized_leak_pot),  "Testing Full net", "0"
        )

def set_init_activity_bioparams(params, base, data_reader="KimDataReaderV3"):
    magnitude, unit = split_neuroml_quantity(base)
    names, conns = get_cell_names_and_connection(data_reader)

    for neuron in names:

        randomized_leak_pot = magnitude + random.normalvariate(0, 1) * 1e-4

        params.add_bioparameter(
            "%s_s_init" % (get_c302_name(neuron)), "%.6f" % (randomized_leak_pot),  "Testing Full net", "0"
        )


def setup(
    parameter_set,
    generate=False,
    duration=10000,
    dt=0.1, #ms
    target_directory="experiments",
    muscles_to_include=[], #None,  # None => All!
    data_reader="KimDataReaderV3", #"KimDataReader", #"KimDataReaderV2", #"KimDataReaderV3", #"VarshneyDataReader", #"SpreadsheetDataReader",
    param_overrides={},
    config_param_overrides={},
    verbose=True,
):
    ParameterisedModel = getattr(
        importlib.import_module("c302.parameters_%s" % parameter_set),
        "ParameterisedModel",
    )
    params = ParameterisedModel()

    # cells_to_stimulate = ["PVDL", "PVDR", "PDEL", "PDER"]
    cells_to_stimulate = ["PLML", "PLMR"]

    cells_to_plot = [
"AVBL", "AVBR", "DB1", "DB2", "DB3", "DB4", "DB5", "DB6", "DB7", "PVCL",
"PVCR", "RIBL", "RIBR", "VB1", "VB2", "VB3", "VB4", "VB5", "VB6", "VB7",
"VB8", "VB9", "VB10", "VB11"
    ]
    # cells_to_plot = ["AVBL"]

    stimuli_magnitude = "5nA"

    params.set_bioparameter(
        "unphysiological_offset_current", stimuli_magnitude, "Testing Full net", "0"

    )
    params.set_bioparameter(
        "unphysiological_offset_current_del", "0 ms", "Testing Full net", "0"
    )
    params.set_bioparameter(
        "unphysiological_offset_current_dur", "10000 ms", "Testing Full net", "0"
    )

    set_vth_bioparams(params, cells_to_stimulate, stimuli_magnitude)

    set_init_potential_bioparams(params, "-35 mV", data_reader=data_reader)

    set_init_activity_bioparams(params, "0", data_reader=data_reader)

    reference = "c302_%s_Full_Interactome" % parameter_set

    cell_names, conns = c302.get_cell_names_and_connection(data_reader)

    nml_doc = None

    # conn_polarity_override = {
    #     r"^DB\d+-DD\d+$": "inh",
    #     r"^VB\d+-VD\d+$": "inh",
    # }

    if generate:
        nml_doc = c302.generate(
            reference,
            params,
            cells_to_plot=cells_to_plot,
            cells_to_stimulate=cells_to_stimulate,
            muscles_to_include=muscles_to_include,
            duration=duration,
            dt=dt,
            # vmin=-72 if parameter_set == "A" else -52,
            # vmax=-48 if parameter_set == "A" else -28,
            target_directory=target_directory,
            param_overrides=param_overrides,
            # conn_polarity_override=conn_polarity_override,
            verbose=verbose,
            data_reader=data_reader,
            # print_connections=True
        )

    return cell_names, cells_to_stimulate, params, muscles_to_include, nml_doc


if __name__ == "__main__":
    # parameter_set = sys.argv[1] if len(sys.argv) == 2 else "A"
    parameter_set="I0"

    setup(parameter_set, generate=True)
