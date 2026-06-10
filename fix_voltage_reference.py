import re
import sys
import os

def patch_neuron_script(input_path, output_path=None):
    if output_path is None:
        output_path = input_path  # patch in place
    
    with open(input_path, 'r') as f:
        content = f.read()
    
    # Fix v recording:
    # m_CELLID_neuron_cell_CELLID[0].v  →  CELLID[0].v(0.5)
    content = re.sub(
        r'm_(\w+)_neuron_cell_\1\[0\]\.v',
        r'\1[0].v(0.5)',
        content
    )

    # Fix v recording for muscle cells:
    # m_GenericMuscleCell_MUSCLEID[0].v  →  MUSCLEID[0].v(0.5)
    content = re.sub(
        r'm_GenericMuscleCell_(\w+)\[0\]\.v',
        r'\1[0].v(0.5)',
        content
    )
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    print(f"Patched: {output_path}")


if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("please enter filepath.")
        exit(-1)

    filepath = sys.argv[1]

    if not os.path.isfile(filepath):
        print("please enter legit filepath.")
        exit(-1)

    output_path = None
    if len(sys.argv) > 2:
        output_path = sys.argv[2]

    patch_neuron_script(filepath, output_path)