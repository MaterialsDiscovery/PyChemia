"""
Created on April 25 2020 
@author: Pedram Tavadze
"""
import os
import numpy as np
from .xml_output import parse_vasprun
from .incar import VaspInput
from ..codes import CodeOutput
from ...core import Structure 
from ...visual import DensityOfStates
from ...crystal.kpoints import KPoints



class VaspXML(CodeOutput):
    
    def __init__(self, filename='vasprun.xml'):
        
        CodeOutput.__init__(self)
        if not os.path.isfile(filename):
            raise ValueError('File not found ' + filename)
        else:
            self.filename = filename
        
        self.spins_dict = {'spin 1':'Spin-up','spin 2':'Spin-down'}
        # self.positions = None
        # self.stress = None
        self.bands = None
        #self.array_sizes = {}
        self.data = self.read()


    def read(self):
        return parse_vasprun(self.filename)
    
    

               
    def _get_dos_total(self):

        spins = list(self.data['general']['dos']['total']['array']['data'].keys())
        energies = np.array(self.data['general']['dos']['total']['array']['data'][spins[0]])[:,0]
        dos_total = {'energies':energies}
        for ispin in spins:
            dos_total[self.spins_dict[ispin]] = np.array(self.data['general']['dos']['total']['array']['data'][ispin])[:,1]
        
        return dos_total,list(dos_total.keys())


    def _get_dos_projected(self,atoms=[]):
        
        if len(atoms) == 0:
            atoms = np.arange(self.initial_structure.natom)
       
        if 'partial' in self.data['general']['dos']:
            dos_projected = {}
            ion_list = ["ion %s"%str(x+1) for x in atoms] # using this name as vasrun.xml uses ion #
            for i in range(len(ion_list)):
                iatom = ion_list[i]
                name = self.initial_structure.symbols[atoms[i]]+str(atoms[i])
                spins = list(self.data['general']['dos']['partial']['array']['data'][iatom].keys())
                energies = np.array(self.data['general']['dos']['partial']['array']['data'][iatom][spins[0]][spins[0]])[:,0]
                dos_projected[name] ={'energies':energies}
                for ispin in spins:
                    dos_projected[name][self.spins_dict[ispin]] =np.array(self.data['general']['dos']['partial']['array']['data'][iatom][ispin][ispin])[:,1:]
            return dos_projected,self.data['general']['dos']['partial']['array']['info']
        else:
            print("This calculation does not include partial density of states")
            return None,None

    
    @property
    def dos_to_dict(self):
        """
        Returns the complete density (total,projected) of states as a python dictionary        
        """
        return {'total':self._get_dos_total(),'projected':self._get_dos_projected()}
    
    @property
    def dos_total(self):
        """
        Returns the total density of states as a pychemia.visual.DensityOfSates object
        """
        dos_total,labels = self._get_dos_total()
        dos_total['energies']-=self.fermi
        return DensityOfStates(np.array([dos_total[x] for x in dos_total]).T, 
                                         title='Total Density Of States',labels=[x.capitalize() for x in labels])

    @property
    def dos_projected(self): 
        """
        Returns the a list of projected density of states as a pychemia.visual.DensityOfSates object
        each element refers to each atom
        """
        ret = []
        atoms = np.arange(self.initial_structure.natom,dtype=int)
        dos_projected,info = self._get_dos_projected(atoms=atoms)
        if dos_projected == None:
            return None
        ndos = len(dos_projected[list(dos_projected.keys())[0]]['energies'])
        norbital = len(info)-1
        nspin = len(dos_projected[list(dos_projected.keys())[0]].keys()) -1
        info[0] = info[0].capitalize()
        labels = []
        labels.append(info[0])
        if nspin >1 :
            for il in info[1:]:
                labels.append(il+'-Up')
            for il in info[1:]:
                labels.append(il+'-Down')
        else : 
            labels = info
        for iatom in dos_projected:
            table = np.zeros(shape = (ndos,norbital*nspin+1))
            table[:,0] = dos_projected[iatom]['energies'] - self.fermi
            start = 1
            for key in dos_projected[iatom]:
                if key == 'energies':
                    continue
                end = start + norbital
                table[:,start:end] = dos_projected[iatom][key]
                start = end 
            temp_dos = DensityOfStates(table,title='Projected Density Of States %s'%iatom,labels=labels)
            ret.append(temp_dos)
        return ret



    def dos_parametric(self,atoms=None,orbitals=None,spin=None,title=None):
        """
        This function sums over the list of atoms and orbitals given 
        for example dos_paramateric(atoms=[0,1,2],orbitals=[1,2,3],spin=[0,1])
        will sum all the projections of atoms 0,1,2 and all the orbitals of 1,2,3 (px,py,pz)
        and return separatly for the 2 spins as a DensityOfStates object from pychemia.visual.DensityofStates
        
        :param atoms: list of atom index needed to be sumed over. count from zero with the same 
                      order as POSCAR
        
        :param orbitals: list of orbitals needed to be sumed over 
        |  s  ||  py ||  pz ||  px || dxy || dyz || dz2 || dxz ||x2-y2||
        |  0  ||  1  ||  2  ||  3  ||  4  ||  5  ||  6  ||  7  ||  8  ||
        
        :param spin: which spins to be included. count from 0
                     There are no sum over spins
        
        """
        projected = self.dos_projected
        if atoms == None :
            atoms = np.arange(self.initial_structure.natom,dtype=int)
        if spin == None :
            spin = [0,1]
        if orbitals == None :
            orbitals = np.arange((len(projected[0].labels)-1)//2,dtype=int)
        if title == None:
            title = 'Sum'
        orbitals = np.array(orbitals)
        if len(spin) ==2:
            labels = ['Energy','Spin-Up','Spin-Down']
            new_orbitals = []
            for ispin in spin :
                new_orbitals.append(list(orbitals+ispin*(len(projected[0].labels)-1)//2))
            orbitals = new_orbitals
        else : 
            if spin[0] == 0:
                labels = ['Energy','Spin-Up']
            elif spin[0] == 1:
                labels = ['Energy','Spin-Down']
            
        
        
        ret = np.zeros(shape=(len(projected[0].energies),len(spin)+1))
        ret[:,0] = projected[0].energies
        
        for iatom in atoms :
            if len(spin) == 2 :
                ret[:,1:]+=self.dos_projected[iatom].values[:,orbitals].sum(axis=2)
            elif len(spin) == 1 :
                ret[:,1]+=self.dos_projected[iatom].values[:,orbitals].sum(axis=1)
                
        return DensityOfStates(table=ret,title=title,labels=labels)
        
    
    def get_band_projection(self):
        return 

    

    @property 
    def kpoints(self):
        """
        Returns the kpoints used in the calculation in form of a pychemia.core.KPoints object
        """
        
        if self.data['kpoints_info']['mode'] == 'listgenerated':
            kpoints = KPoints(kmode='path',kvertices=self.data['kpoints_info']['kpoint_vertices'])
        else :
            kpoints = KPoints(kmode=self.data['kpoints_info']['mode'].lower(),
                               grid=self.data['kpoints_info']['kgrid'],
                               shifts=self.data['kpoints_info']['user_shift'])
        return kpoints
    
    @property
    def kpoints_list(self):
        """
        Returns the list of kpoints and weights used in the calculation in form of a pychemia.core.KPoints object
        """
        return KPoints(kmode='reduced',kpoints_list=self.data['kpoints']['kpoints_list'], 
                       weights=self.data['kpoints']['k_weights'])  
    
    @property 
    def incar(self):
        """
        Returns the incar parameters used in the calculation as pychemia.code.vasp.VaspIncar object
        """
        return VaspInput(variables=self.data['incar'])
        
    @property
    def final_data(self):
        """
        Returns the final free energy, energy_wo_entropy and energy_sigma>0 as a python dictionary
        """
        return {'energy':{'free_energy':self.iteration_data[-1]['energy']['e_fr_energy'],
                          'energy_without_entropy':self.iteration_data[-1]['energy']['e_wo_entrp'],
                          'energy(sigma->0)':self.iteration_data[-1]['energy']['e_0_energy']}}

    @property 
    def vasp_parameters(self):
        """
        Returns all of the parameters vasp has used in this calculation
        """
        return self.data['vasp_params']

    @property 
    def potcar_info(self):
        """
        Returns the information about pseudopotentials(POTCAR) used in this calculation
        """
        return self.data['atom_info']['atom_types']

    @property 
    def fermi(self):
        """
        Returns the fermi energy
        """
        return self.data['general']['dos']['efermi']

    @property
    def species(self):
        """
        Returns the species in POSCAR
        """
        return self.initial_structure.species

    @property 
    def structures(self):
        """
        Returns a list of pychemia.core.Structure representing all the ionic step structures
        """
        symbols = [x.strip() for x in self.data['atom_info']['symbols']]
        structures = []
        for ist in self.data['structures']:
            structures.append(Structure(symbols=symbols,reduced=ist['reduced'],cell=ist['cell']))
        return structures
    
    
    @property 
    def forces(self):
        """
        Returns all the forces in ionic steps
        """
        return self.data['forces']
    
    
    @property
    def initial_structure(self):
        """
        Returns the initial Structure as a pychemia structure
        """
        return self.structures[0]

    @property 
    def final_structure(self):
        """
        Returns the final Structure as a pychemia structure
        """
        
        return self.structures[-1]

    @property 
    def iteration_data(self):
        """
        Returns a list of information in each electronic and ionic step of calculation
        """
        return self.data['calculation']
    
    @property 
    def energies(self):
        """
        Returns a list of energies in each electronic and ionic step [ionic step,electronic step, energy]
        """
        scf_step = 0
        ion_step = 0
        double_counter = 1
        energies = []
        for calc in self.data['calculation']:
            if 'ewald' in calc['energy']:
                if double_counter == 0 :
                    double_counter+=1
                    scf_step +=1
                elif double_counter == 1:
                    double_counter = 0
                    ion_step += 1
                    scf_step = 1  
            else :
                scf_step += 1
            energies.append([ion_step,scf_step,calc['energy']['e_0_energy']])
        return energies

    @property
    def last_energy(self):
        """
        Returns the last calculated energy of the system
        """
        return self.energies[-1][-1]

    @property 
    def energy(self):
        """
        Returns the last calculated energy of the system
        """
        return self.last_energy

    @property
    def convergence_electronic(self):
        """
        Returns a boolian representing if the last electronic self-consistent 
        calculation converged
        """
        ediff = self.vasp_parameters['electronic']['EDIFF']
        last_dE = abs(self.energies[-1][-1]-self.energies[-2][-1])
        if last_dE < ediff :
            return True
        else :
            return False
        
    @property
    def convergence_ionic(self):
        """
        Returns a boolian representing if the ionic part of the  
        calculation converged
        """
        energies = np.array(self.energies)
        nsteps = len(np.unique(np.array(self.energies)[:,0]))
        if nsteps == 1:
            print('This calculation does not have ionic steps')
            return True
        else :
            ediffg = self.vasp_parameters['ionic']['EDIFFG']
            if ediffg < 0 :
                last_forces_abs = np.abs(np.array(self.forces[-1]))
                return not(np.any(last_forces_abs > abs(ediffg)))
            else :
                last_ionic_energy = energies[(energies[:,0] == nsteps)][-1][-1]
                penultimate_ionic_energy = energies[(energies[:,0] == (nsteps-1))][-1][-1]
                last_dE = abs(last_ionic_energy-penultimate_ionic_energy)
                if last_dE < ediffg:
                    return True
        return False
            
    @property
    def convergence(self):
        """
        Returns a boolian representing if the the electronic self-consistent 
        and ionic calculation converged
        """
        return (self.convergence_electronic and self.convergence_ionic)
    
    
    @property
    def is_finished(self):
        """
        Always returns True, need to fix this according to reading the xml as if the calc is 
        not finished we will have errors in xml parser
        """
        # if vasprun.xml is read the calculation is finished
        return True