import os
import shutil
import json
import logging
from pychemia.code.vasp import VaspJob, VaspOutput, VaspInput
from pychemia.crystal import KPoints
from pychemia.code.vasp import read_poscar
from pychemia.utils.mathematics import round_small

__author__ = 'Guillermo Avendano-Franco'


class RelaxPopulation:
    def __init__(self, population, basedir, target_force=1E-2, target_stress=1E-2):
        self.population = population
        self.basedir = basedir
        self.vasp_jobs = {}
        self.runs = {}
        self.runner = None
        self.status = {}
        self.target_force = target_force
        self.target_stress = target_stress

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.basedir):
            os.makedirs(self.basedir)
        elif clean:
            for i in os.listdir(self.basedir):
                shutil.rmtree(self.basedir + os.sep + i)
        for i in self.population.pcdb.entries.find():
            name = self.basedir + os.sep + str(i['_id'])
            if not os.path.isdir(name):
                os.mkdir(name)

    def create_inputs(self, kp_density=10000, encut=1.0):
        # kpoints = KPoints(kmode='gamma', grid=[4, 4, 4])
        for entry in self.population.pcdb.entries.find():
            name = str(entry['_id'])
            workdir = self.basedir + os.sep + name
            structure = self.population.db.get_structure(entry['_id'])
            kpoints = KPoints.optimized_grid(structure.lattice, kp_density=kp_density)
            print(kpoints)
            vj = VaspJob(workdir=workdir)
            vj.initialize(structure=structure, kpoints=kpoints)
            inp = VaspInput()
            inp.set_rough_relaxation()
            vj.set_input_variables(inp)
            vj.write_potcar()
            vj.input_variables.set_encut(ENCUT=encut, POTCAR=workdir + os.sep + 'POTCAR')
            vj.set_inputs()
            self.vasp_jobs[name] = vj
            self.runs[name] = 0
            self.status[name] = ['ACTIVE']

    def add_status(self, entry_id, value):
        if value not in self.status[entry_id]:
            self.status[entry_id].append(value)

    def del_status(self, entry_id, value):
        if value in self.status[entry_id]:
            self.status[entry_id].remove(value)

    def flip_status(self, entry_id, oldvalue, newvalue):
        self.del_status(entry_id, oldvalue)
        self.add_status(entry_id, newvalue)

    def modify_input(self, entry_id):
        if 'RELAXED' not in self.status[entry_id] and 'NOPROCAR' not in self.status[entry_id] \
                and 'NOOUTCAR' not in self.status[entry_id]:
            return True
        else:
            return False

    def update(self, workdir):
        """
        This routine determines how to proceed with the relaxation
        for one specific work directory

        :param workdir: (str) String representation of the id in the mongodb
        :return:
        """

        # workdir = self.basedir + os.sep + entry_id
        entry_id = os.path.basename(workdir)
        vj = self.vasp_jobs[entry_id]
        runj = self.runs[entry_id]
        if os.path.isfile(workdir + os.sep + 'OUTCAR'):
            vj.get_outputs()
        self.update_history(entry_id)

        if os.path.isfile(workdir + os.sep + 'RELAXED'):
            self.add_status(entry_id, 'RELAXED')
        elif not os.path.isfile(workdir + os.sep + 'PROCAR'):
            self.add_status(entry_id, 'NOPROCAR')
        else:
            self.del_status(entry_id, 'NOPROCAR')
            if not os.path.isfile(workdir + os.sep + 'OUTCAR'):
                self.add_status(entry_id, 'NOOUTCAR')
            else:
                self.del_status(entry_id, 'NOOUTCAR')
                print('-')
                vo = VaspOutput(workdir + os.sep + 'OUTCAR')
                relaxation_info = vo.relaxation_info()
                if len(relaxation_info) != 3:
                    print('[' + str(entry_id) + ']' + ' Missing some data in OUTCAR (forces or stress)')
                    self.add_status(entry_id, 'NOOUTCAR')

                print('[' + str(entry_id) + ']' + 'Results:')
                for i in relaxation_info:
                    print('[' + str(entry_id) + '] %20s %12.5e' % (i, relaxation_info[i]))

                # Conditions to consider the structure relaxed
                if relaxation_info['avg_force'] < self.target_force:
                    if relaxation_info['avg_stress_diag'] < self.target_stress:
                        if relaxation_info['avg_stress_non_diag'] < self.target_stress:
                            wf = open(workdir + os.sep + 'RELAXED', 'w')
                            for i in relaxation_info:
                                wf.write("%15s %12.3f" % (i, relaxation_info[i]))
                            wf.close()
                            wf = open(workdir + os.sep + 'COMPLETE', 'w')
                            for i in relaxation_info:
                                wf.write("%15s %12.3f" % (i, relaxation_info[i]))
                            wf.close()
                            self.add_status(entry_id, 'RELAXED')

        if self.modify_input(entry_id):

            # How to change ISIF
            if relaxation_info['avg_force'] < 0.1:
                if relaxation_info['avg_stress_diag'] < 0.1:
                    if relaxation_info['avg_stress_non_diag'] < 0.1:
                        vj.input_variables.variables['ISIF'] = 3
                    else:
                        vj.input_variables.variables['ISIF'] = 3
                else:
                    vj.input_variables.variables['ISIF'] = 3
            else:
                vj.input_variables.variables['ISIF'] = 2

            # How to change IBRION
            # if info['avg_force'] < 0.1 and info['avg_stress_diag'] < 0.1 and info['avg_stress_non_diag'] < 0.1:
            #    vj.input_variables.variables['IBRION'] = 1
            # elif info['avg_force'] < 1 and info['avg_stress_diag'] < 1 and info['avg_stress_non_diag'] < 1:
            #    vj.input_variables.variables['IBRION'] = 2
            # else:
            #    vj.input_variables.variables['IBRION'] = 3

            # How to change EDIFF
            if vj.input_variables.variables['EDIFF'] > 2 * 1E-4:
                vj.input_variables.variables['EDIFF'] = round_small(vj.input_variables.variables['EDIFF'] / 2)
            else:
                vj.input_variables.variables['EDIFF'] = 1E-4

            # How to change EDIFFG
            if vj.input_variables.variables['EDIFFG'] < - 2 * self.target_force:
                vj.input_variables.variables['EDIFFG'] = round_small(vj.input_variables.variables['EDIFFG'] / 2)
            else:
                vj.input_variables.variables['EDIFFG'] = - self.target_force

            # Print new values
            print('[' + str(entry_id) + ']' + 'New Values:')
            for i in ['ISIF', 'IBRION', 'EDIFF', 'EDIFFG']:
                print('[' + str(entry_id) + ']' + i + ' : ', vj.input_variables.variables[i])
            print('-')

            for i in ['OUTCAR']:
                if not os.path.exists(workdir + os.sep + i):
                    wf = open(workdir + os.sep + i, 'w')
                    wf.write('')
                    wf.close()
                log = logging.handlers.RotatingFileHandler(workdir + os.sep + i, maxBytes=1, backupCount=1000)
                log.doRollover()

            try:
                vj.structure = read_poscar(workdir + os.sep + 'CONTCAR')
            except ValueError:
                print('Error reading CONTCAR')

            vj.set_inputs()
            properties = vj.outcar
            status = self.status[entry_id]
            newentry = self.population.db.update(entry_id, structure=vj.structure, properties=properties, status=status)

            vj.save_json(workdir + os.sep + 'vaspjob.json')
            wf = open(workdir + os.sep + 'entry.json', 'w')
            json.dump(newentry, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
            return True
        else:
            vj.set_inputs()
            status = self.status[entry_id]
            newentry = self.population.db.update(entry_id, structure=vj.structure, status=status)

            vj.save_json(workdir + os.sep + 'vaspjob.json')
            wf = open(workdir + os.sep + 'entry.json', 'w')
            json.dump(newentry, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
            return True

    def update_history(self, entry_id):
        filename = 'pychemia_relaxation.json'
        filepath = self.basedir + os.sep + entry_id + os.sep + filename
        if not os.path.exists(filepath):
            wf = open(filepath, 'w')
            data = [self.vasp_jobs[entry_id].to_dict]
            json.dump(data, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
        else:
            rf = open(filepath, 'r')
            data = json.load(rf)
            rf.close()
            data.append(self.vasp_jobs[entry_id].to_dict)
            wf = open(filepath, 'w')
            json.dump(data, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()

    @property
    def workdirs(self):
        return [self.basedir + os.sep + name for name in self.population.members]

    @property
    def active_workdirs(self):
        return [self.basedir + os.sep + name for name in self.population.actives]

    def run(self, runner):

        entries_ids = self.population.members

        def worker(workdir):
            wf = open(workdir + os.sep + 'LOCK', 'w')
            wf.write('')
            wf.close()
            runner.run()
            os.remove(workdir + os.sep + 'LOCK')

        def checker(workdir):
            if os.path.isfile(workdir + os.sep + 'LOCK'):
                return False
            return self.update(workdir)

        workdirs = [self.basedir + os.sep + i for i in self.population.actives]
        runner.run_multidirs(workdirs, worker, checker)

        if not self.is_running:
            self.run(runner)

    def set_run(self, code, runner, basedir, kp_density=10000, encut=1.1):

        self.runner = runner

        self.create_dirs(clean=True)
        self.create_inputs(kp_density=kp_density, encut=encut)
