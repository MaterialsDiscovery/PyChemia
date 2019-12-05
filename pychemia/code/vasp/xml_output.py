"""
Created on Mon Dec 02 2019 
@author: Pedram Tavadze
"""
import pychemia                     
import xml.etree.ElementTree as ET    
import numpy as np # we ddon't need this as well

def TextToBool(text) :
    '''boolians in vaspxml are stores as T or F in str format, this function coverts them to python boolians '''
    text = text.strip(' ')
    if text == 'T' or text == '.True.' or text == '.TRUE.':
        return True
    else :
        return False
    
def conv(ele,_type): 
    '''This function converts the xml text to the type specified in the attrib of xml tree '''
    if   _type == 'string' : 
        return str(ele) 
    elif _type == 'int' : 
        return int(ele) 
    elif _type == 'logical': 
        return TextToBool(ele) 
    elif _type == 'float' : 
        return float(ele) 

def get_varray(xml_tree) : 
    '''Returns an array for each varray tag in vaspxml '''
    ret = [] 
    for ielement in xml_tree: 
        ret.append([float(x) for x in ielement.text.split()]) 
    return ret

def get_params(xml_tree,dest): 
    '''dest should be a dictionary 
    This function is recurcive #check spelling''' 
    for ielement in xml_tree : 
        if ielement.tag =='separator' : 
            dest[ielement.attrib['name']] = {} 
            dest[ielement.attrib['name']] = get_params(ielement, dest[ielement.attrib['name']]) 
        else : 
            if 'type' in ielement.attrib : 
                _type =  ielement.attrib['type'] 
            else : 
                _type = 'float' 
            if ielement.text == None : 
                dest[ielement.attrib['name']] = None 
            elif len(ielement.text.split()) > 1 : 
                dest[ielement.attrib['name']] = [ conv(x,_type) for x in ielement.text.split() ] 
            else : 
                dest[ielement.attrib['name']] = conv(ielement.text,_type) 
    return dest 

def get_structure(xml_tree):
    '''Returns a dictionary of the structure '''
    ret = {}
    for ielement in xml_tree : 
        if ielement.tag == 'crystal': 
            for isub in ielement: 
                if isub.attrib['name'] == 'basis': 
                    ret['cell'] = get_varray(isub)
                elif isub.attrib['name'] == 'volume': 
                    ret['volume'] = float(isub.text) 
                elif isub.attrib['name'] == 'rec_basis': 
                    ret['rec_cell'] = get_varray(isub)
        elif ielement.tag == 'varray' : 
            if ielement.attrib['name'] == 'positions': 
                ret['reduced'] = get_varray(ielement)
    return ret

def get_scstep(xml_tree) :
    '''This function extracts the self-consistent step information '''
    scstep = {} 
    scstep['time'] = {} # Do we need this? 
    scstep['energy'] = {} 
    for isub in xml_tree : 
        if isub.tag=='time': 
            scstep['time'][isub.attrib['name']] = [float(x) for x in isub.text.split()] 
        elif isub.tag == 'energy': 
            for ienergy in isub : 
                scstep['energy'][ienergy.attrib['name']] = float(ienergy.text) 
    return scstep



def get_set(xml_tree,ret): 
    ''' This function will extract any element taged set recurcively'''
    if xml_tree[0].tag == 'r': 
        ret[xml_tree.attrib['comment']] = get_varray(xml_tree)         
        return ret 
    else :     
        for ielement in xml_tree: 
            if ielement.tag == 'set': 
                ret[ielement.attrib['comment']] = {} 
                ret[ielement.attrib['comment']] = get_set(ielement,ret[ielement.attrib['comment']])
        return ret 
        
def get_general(xml_tree,ret): 
    ''' This function will parse any element in calculatio other than the structures, scsteps''' 
    if 'dimension' in [x.tag for x in xml_tree]: 
        ret['info'] = [] 
        ret['data'] = {} 
        for ielement in xml_tree: 
            if ielement.tag == 'field': 
                ret['info'].append(ielement.text.strip(' ')) 
            elif ielement.tag == 'set': 
                for iset in ielement: 
                    ret['data'] = get_set(iset,ret['data']) 
        return ret 
    else : 
        for ielement in xml_tree: 
            if ielement.tag == 'i': 
                if 'name' in ielement.attrib :
                    if ielement.attrib['name']=='efermi':
                        ret['efermi'] = float(ielement.text)
                continue 
            ret[ielement.tag] = {} 
            ret[ielement.tag] = get_general(ielement, ret[ielement.tag]) 
        return ret   


def parse_vasprun(vasprun):
    tree = ET.parse(vasprun)
    root = tree.getroot()

    calculation = []
    structures  = []
    forces      = []
    stresses    = []
    orbital_magnetization = {}
    run_info = {}
    incar = {}
    general = {}
    kpoints_info = {}
    vasp_params = {}
    kpoints_list = [] 
    k_weights = []
    atom_info = {} 
    for ichild in root:
    
        
        if ichild.tag == 'generator':
            for ielement in ichild:
                run_info[ielement.attrib['name']] = ielement.text
                
                
                
        elif ichild.tag == 'incar':
            incar = get_params(ichild,incar)
    
        ## Skipping 1st structure which is primitive cell
    
        elif ichild.tag == 'kpoints' : 
    
            for ielement in ichild :
                if ielement.items()[0][0] == 'param': 
                    kpoints_info['mode'] = ielement.items()[0][1] 
                    for isub in ielement :  
                        if isub.attrib['name'] == 'divisions': 
                            kpoints_info['kgrid'] = [int(x) for x in isub.text.split()] 
                        elif isub.attrib['name'] == 'usershift': 
                            kpoints_info['user_shift'] = [float(x) for x in isub.text.split()] 
                        elif isub.attrib['name'] == 'genvec1': 
                            kpoints_info['genvec1'] = [float(x) for x in isub.text.split()] 
                        elif isub.attrib['name'] == 'genvec2': 
                            kpoints_info['genvec2'] = [float(x) for x in isub.text.split()] 
                        elif isub.attrib['name'] == 'genvec3': 
                            kpoints_info['genvec3'] = [float(x) for x in isub.text.split()] 
                        elif isub.attrib['name'] == 'shift': 
                            kpoints_info['shift'] = [float(x) for x in isub.text.split()] 
    
                elif ielement.items()[0][1] == 'kpointlist': 
                    for ik in ielement : 
                        kpoints_list.append([float(x) for x in ik.text.split()]) 
                    kpoints_list = np.array(kpoints_list)
                elif ielement.items()[0][1] == 'weights': 
                    for ik in ielement :
                        k_weights.append(float(ik.text))
                    k_weights = np.array(k_weights)
        
    
                    
        ## Vasp Parameters 
        elif ichild.tag == 'parameters' : 

            vasp_params = get_params(ichild,vasp_params) 
    
        ## Atom info
    
        elif ichild.tag == 'atominfo' : 

            for ielement in ichild :  
                if ielement.tag == 'atoms' : 
                    atom_info['natom'] = int(ielement.text) 
                elif ielement.tag == 'types' : 
                    atom_info['nspecies'] = int(ielement.text) 
                elif ielement.tag == 'array' : 
                    if ielement.attrib['name'] == 'atoms': 
                        for isub in ielement : 
                            if isub.tag == 'set': 
                                atom_info['symbols'] = [] 
                                for isym  in isub : 
                                    atom_info['symbols'].append(isym[0].text) 
                    elif ielement.attrib['name'] == 'atomtypes': 
                        atom_info['atom_types'] = {} 
                        for isub in ielement : 
                            if isub.tag == 'set': 
                                for iatom in isub : 
                                    atom_info['atom_types'][iatom[1].text] = {} 
                                    atom_info['atom_types'][iatom[1].text]['natom_per_specie'] = int(iatom[0].text) 
                                    atom_info['atom_types'][iatom[1].text]['mass'] = float(iatom[2].text) 
                                    atom_info['atom_types'][iatom[1].text]['valance'] = float(iatom[3].text) 
                                    atom_info['atom_types'][iatom[1].text]['pseudopotential'] = iatom[4].text 
    
    
    
        elif ichild.tag == 'structure' : 
            if ichild.attrib['name'] == 'initialpos': 
                initial_pos = get_structure(ichild)
            elif ichild.attrib['name'] == 'finalpos': 
                final_pos = get_structure(ichild)
                
                
                
        
        
        elif ichild.tag == 'calculation' : 
            for ielement in ichild : 
                if  ielement.tag == 'scstep': 
                    calculation.append(get_scstep(ielement)) 
                elif ielement.tag == 'structure' :  
                    structures.append(get_structure(ielement))
                elif ielement.tag == 'varray': 
                    if ielement.attrib['name'] == 'forces': 
                        forces.append(get_varray(ielement)) 
                    elif ielement.attrib['name'] == 'stress': 
                        stresses.append(get_varray(ielement)) 
                        
                # elif ielement.tag == 'eigenvalues': 
                #     for isub in ielement[0] :  
                #         if isub.tag == 'set': 
                #             for iset in isub : 
                #                 eigen_values[iset.attrib['comment']] = {} 
                #                 for ikpt in iset :  
                #                     eigen_values[iset.attrib['comment']][ikpt.attrib['comment']] = get_varray(ikpt) 
                
                elif ielement.tag == 'separator': 
                    if ielement.attrib['name'] == "orbital magnetization": 
                        for isub in ielement :  
                            orbital_magnetization[isub.attrib['name']] = [float(x) for x in isub.text.split()]
                            
                            
                # elif ielement.tag == 'dos': 
                #     for isub in ielement :  
                #         if 'name' in isub.attrib: 
                #             if isub.attrib['name'] == 'efermi' :  
                #                 dos['efermi'] = float(isub.text) 
                #             else :  
                #                 dos[isub.tag] = {} 
                #                 dos[isub.tag]['info'] = [] 
                #               for iset in isub[0]  : 
                #                   if iset.tag == 'set' : 
                #                       for isub_set in iset: 
                #                           dos[isub.tag] = get_set(isub_set,dos[isub.tag]) 
                #                   elif iset.tag == 'field' : 
                #                       dos[isub.tag]['info'].append(iset.text.strip(' '))
                else :
                    general[ielement.tag] = {}
                    general[ielement.tag] = get_general(ielement,general[ielement.tag])
                    
                    
    return {'calculation':calculation,'structures':structures,'forces':forces,'run_info':run_info,'incar':incar,'general':general,'kpoints_info':kpoints_info,'vasp_params':vasp_params,'kpoints':{'kpoints_list':kpoints_list,'k_weights':k_weights},'atom_info':atom_info}


if __name__ == '__main__':
    vi = pychemia.code.vasp.VaspInput(variables=incar)   
    kp_fromGen = pychemia.crystal.KPoints(kmode=kpoints_info['mode'].lower(),shifts=kpoints_info['user_shift'],grid=kpoints_info['kgrid'])
    kp_fromList = pychemia.crystal.KPoints(kpoints_list=kpoints_list,weights=k_weights) # Ask Guillermo why this does not give the correct answer
    st_init = pychemia.core.Structure(symbols=atom_info['symbols'],reduced=initial_pos['reduced'],cell=initial_pos['cell'])
    st_final = pychemia.core.Structure(symbols=atom_info['symbols'],reduced=final_pos['reduced'],cell=final_pos['cell'])
