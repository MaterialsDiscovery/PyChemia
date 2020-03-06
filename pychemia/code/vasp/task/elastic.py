import os
import re
import time
import json
import numpy as np
from pychemia.crystal import KPoints
from pychemia import pcm_log
from pychemia.utils.serializer import generic_serializer
from ..vasp import VaspJob
from ..outcar import read_vasp_stdout
from ...tasks import Task

__author__ = 'Guillermo Avendano-Franco'


class ElasticModuli(Task):
    def __init__(self, structure, workdir='.', executable='vasp', encut=1.3, kpoints=None, kp_density=1E4):

        self.encut = encut
        if kpoints is None:
            kp = KPoints.optimized_grid(structure.lattice, kp_density=kp_density, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kpoints
        self.task_params = {'encut': self.encut, 'kpoints': self.kpoints.to_dict}
        Task.__init__(self, structure=structure, task_params=self.task_params, workdir=workdir, executable=executable)

    def run(self, nparal=4):

        vj = VaspJob(workdir=self.workdir, executable=self.executable)
        vj.initialize(self.structure, self.kpoints)
        vj.clean()
        vj.job_static()
        vj.input_variables.set_density_for_restart()
        vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.input_variables.variables['IBRION'] = 6
        vj.input_variables.variables['ISIF'] = 3
        vj.set_inputs()
        self.encut = vj.input_variables.variables['ENCUT']
        vj.run(mpi_num_procs=nparal)
        pcm_log.debug('Starting VASP')
        while True:
            energy_str = ''
            filename = self.workdir + os.sep + 'vasp_stdout.log'
            if os.path.exists(filename):
                vasp_stdout = read_vasp_stdout(filename=filename)
                if len(vasp_stdout['data']) > 2:
                    scf_energies = [i[2] for i in vasp_stdout['data']]
                    energy_str = ' %7.3f' % scf_energies[1]
                    for i in range(1, len(scf_energies)):
                        if scf_energies[i] < scf_energies[i - 1]:
                            energy_str += ' >'
                        else:
                            energy_str += ' <'
                    pcm_log.debug(energy_str)

            if vj.runner is not None and vj.runner.poll() is not None:
                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    vasp_stdout = read_vasp_stdout(filename=filename)
                    if len(vasp_stdout['data']) > 2:
                        scf_energies = [i[2] for i in vasp_stdout['data']]
                        energy_str += ' %7.3f' % scf_energies[-1]
                        pcm_log.debug(energy_str)
                pcm_log.debug('Execution complete')
                break
            time.sleep(5)
        vj.get_outputs()

        self.output = self.get_elastic_moduli()
        if vj.outcar.is_finished:
            self.finished = True

    def plot(self, filedir=None, file_format='pdf'):
        if filedir is None:
            filedir = self.workdir
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5))
        plt.subplots_adjust(left=0.04, bottom=0.03, right=0.98, top=0.98, wspace=0.1, hspace=0.1)

        data = np.array(self.output['total_elastic_moduli'])
        ax1.imshow(data, interpolation='none', cmap=plt.get_cmap('Greys'))

        for i in range(len(data)):
            for j in range(len(data[i])):
                if data[i, j] > 1E-10:
                    if data[i, j] > 0.75 * np.max(data.flatten()):
                        color = 'white'
                    else:
                        color = 'black'
                    ax1.text(i, j, int(data[i, j]), ha='center', va='center', color=color)

        ax1.set_xticks(np.arange(6))
        ax1.set_xticklabels(['XX', 'YY', 'ZZ', 'XY', 'YZ', 'ZX'])
        ax1.set_yticks(np.arange(6))
        ax1.set_yticklabels(['XX', 'YY', 'ZZ', 'XY', 'YZ', 'ZX'])
        ax1.set_xlabel('Total Elastic Moduli')
        ax1.xaxis.tick_top()

        data = np.array(self.output['elastic_moduli_contributions_from_ionic_relaxation'])
        maxdata = np.abs(np.max(data.flatten()))
        ax2.imshow(data, interpolation='none', cmap=plt.get_cmap('seismic'), vmin=-maxdata, vmax=maxdata)

        for i in range(len(data)):
            for j in range(len(data[i])):
                if np.abs(data[i, j]) > 1E-10:
                    if np.abs(data[i, j]) > 0.75 * maxdata:
                        color = 'white'
                    else:
                        color = 'black'
                    ax2.text(i, j, int(data[i, j]), ha='center', va='center', color=color)

        ax2.set_xticks(np.arange(6))
        ax2.set_xticklabels(['XX', 'YY', 'ZZ', 'XY', 'YZ', 'ZX'])
        ax2.set_yticks(np.arange(6))
        ax2.set_yticklabels(['XX', 'YY', 'ZZ', 'XY', 'YZ', 'ZX'])
        ax2.set_xlabel('Elastic Moduli contributions from Ion relax')
        ax2.xaxis.tick_top()

        plt.savefig(filedir + os.sep + 'elastic.' + file_format)
        return plt.gcf()

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        rf = open(filename)
        data = json.load(rf)
        rf.close()
        self.task_params = data['task_params']
        self.output = data['output']
        self.encut = self.task_params['encut']
        self.kpoints = KPoints.from_dict(self.task_params['kpoints'])

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E
        self.plot(filedir=self.workdir + os.sep + 'REPORT', file_format='jpg')

        style = """
table {
    width:70%;
}
table, th, td {
    border: 1px solid black;
    border-collapse: collapse;
}
th, td {
    padding: 5px;
    text-align: left;
}
table#t01 tr:nth-child(even) {
    background-color: #eee;
}
table#t01 tr:nth-child(odd) {
   background-color:#fff;
}
table#t01 th	{
    background-color: black;
    color: white;
}
"""
        mech = self.get_mechanical_properties()
        table = []
        for i in mech:
            table.append(E.tr(E.td(i), E.td(mech[i]['units']),
                              E.td('%7.3f' % mech[i]['Voigt']),
                              E.td('%7.3f' % mech[i]['Reuss']),
                              E.td('%7.3f' % (0.5 * (mech[i]['Voigt'] + mech[i]['Reuss'])))))

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP Elastic Moduli"), E.style(style)),
                                  E.body(E.h1("VASP Elastic Moduli"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Elastic Moduli'),
                                         E.p(E.img(src='elastic.jpg', width="900", height="500",
                                                   alt="Static Calculation")),
                                         E.table(E.tr(E.th('Property'), E.th('Units'), E.th('Voigt'), E.th('Reuss'),
                                                      E.th('Average')), E.tr(*tuple(table)), id="t01"), E.br()
                                         ))

        return self.report_end(html, file_format)

    def get_elastic_moduli(self):

        if not os.path.isfile(self.workdir + os.sep + 'OUTCAR'):
            print('OUTCAR file not found')
            return

        return elastic_moduli()

    def get_mechanical_properties(self):
        """
        Return the mechanical properties as a dictionary (Total Elastic Moduli)

        :rtype : dict
        """
        return mechanical_properties(self.get_elastic_moduli()['total_elastic_moduli'])


def elastic_moduli(filename='OUTCAR'):
    ret = {}
    rf = open(filename)
    data = rf.read()

    subblock = re.findall('ELASTIC MODULI CONTR FROM IONIC RELAXATION \(kBar\)[\s\d\w]*-*([\s\d\w.-]*)\n\n\n',
                          data)
    ret['elastic_moduli_contr'] = np.array(np.array(subblock[0].split()[:42]).reshape(6, -1)[:, 1:], dtype=float)

    subblock = re.findall('TOTAL ELASTIC MODULI \(kBar\)[\s\d\w]*-*([\s\d\w.-]*)\n\n\n', data)
    ret['total_elastic_moduli'] = np.array(np.array(subblock[0].split()[:42]).reshape(6, -1)[:, 1:], dtype=float)

    return generic_serializer(ret)


def mechanical_properties(ls_elastic_moduli):
    np_elastic_moduli = np.array(ls_elastic_moduli)

    # Forcing the matrix to be symmetric
    np_elastic_moduli = 0.5 * (np_elastic_moduli + np.transpose(np_elastic_moduli))

    elastic_moduli_inv = np.linalg.inv(np_elastic_moduli)

    # Voigt
    em = np_elastic_moduli[:3, :3]
    kkv = np.trace(em) / 90.0 + 2 * np.sum(np.triu(em, 1)) / 90.0
    ggv = np.trace(em) - np.sum(np.triu(em, 1)) + 3 * np.trace(np_elastic_moduli[3:6, 3:6])
    ggv /= 150

    eev = 1.0 / 3.0 / ggv + 1.0 / 9.0 / kkv
    eev = 1.0 / eev

    vv = 3.0 * ggv / (3.0 * kkv + ggv)
    vv = 1.0 - vv
    vv /= 2

    # Reuss
    kkr = np.trace(elastic_moduli_inv[:3, :3]) + 2 * np.sum(np.triu(elastic_moduli_inv[:3, :3], 1))
    kkr = 0.1 / kkr

    eminv = elastic_moduli_inv[:3, :3]
    ggr = 4 * np.trace(eminv) - 4 * np.sum(np.triu(eminv, 1)) + 3 * np.trace(elastic_moduli_inv[3:6, 3:6])
    ggr = 1.5 / ggr

    eer = 1.0 / 3.0 / ggr + 1.0 / 9.0 / kkr
    eer = 1.0 / eer

    vr = 3.0 * ggr / (3.0 * kkr + ggr)
    vr = 1.0 - vr
    vr /= 2

    ret = {'Bulk modulus': {'units': 'GPa', 'Voigt': kkv, 'Reuss': kkr},
           'Shear modulus': {'units': 'GPa', 'Voigt': ggv, 'Reuss': ggr},
           'Young modulus': {'units': 'GPa', 'Voigt': eev, 'Reuss': eer},
           'Poisson ratio': {'units': 'GPa', 'Voigt': vv, 'Reuss': vr}}

    return ret
