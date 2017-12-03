import os
from .input import VaspInput
from .poscar import write_poscar


def write_from_queue(queue, entry_id, destination=None):
    if destination is None:
        dest = '.'
    elif os.path.isfile(destination):
        dest = os.path.dirname(os.path.abspath(destination))
    elif not os.path.exists(destination):
        os.mkdir(destination)
        dest = destination
    elif os.path.isdir(destination):
        dest = destination
    else:
        raise ValueError('Destination not valid')

    st = queue.get_input_structure(entry_id)
    inpvars = queue.get_input_variables(entry_id)

    vi = VaspInput()
    for i in inpvars:
        vi[i] = inpvars[i]

    write_poscar(st, dest + os.sep + 'POSCAR')
    vi.write(dest + os.sep + 'INCAR')

    queue.write_input_files(entry_id, dest)
