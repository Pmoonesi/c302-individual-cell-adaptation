import c302
import sys

import importlib


def setup(
    parameter_set,
    generate=False,
    duration=500,
    dt=0.05, #ms
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

#     cells_to_plot = [
# "AVBL", "AVBR", "DB1", "DB2", "DB3", "DB4", "DB5", "DB6", "DB7", "PVCL",
# "PVCR", "RIBL", "RIBR", "VB1", "VB2", "VB3", "VB4", "VB5", "VB6", "VB7",
# "VB8", "VB9", "VB10", "VB11"
#     ]
    cells_to_plot = ["AVBL"]



    params.set_bioparameter(
        "unphysiological_offset_current", "5pA", "Testing Full net", "0"

    )
    params.set_bioparameter(
        "unphysiological_offset_current_del", "0 ms", "Testing Full net", "0"
    )
    params.set_bioparameter(
        "unphysiological_offset_current_dur", "400 ms", "Testing Full net", "0"
    )

    params.add_bioparameter(
        "AVBL_v_threshold", "10 mV",  "Testing Full net", "0"
    )

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
