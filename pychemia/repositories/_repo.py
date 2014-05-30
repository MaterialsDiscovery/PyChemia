"""
There are two kinds of Repositories in PyChemia

Structure Repositories where many structures are stored
Execution Repositories where  the out of every calculation
is stored

Each structure contains some metadata that is accessible with
the StructureEntry object
Also each calculation has it own metadata accessible by ExecutionEntry
object
"""

import json as _json
import os as _os
import uuid as _uuid
import shutil as _shutil
from pychemia.geometry import load_structure_json
from pychemia.utils.computing import unicode2string


class StructureEntry():
    """
    Defines one entry in the repository of Structures
    """

    def __init__(self, structure=None, repository=None, identifier=None, original_file=None, tags=None):
        """
        Creates a new Entry for Structures
        If identifier is provided the corresponding Structure is load in the Entry
        Otherwise a new entry is created with a UUID random identifier

        Args:
        identifier: (string) UUID identifier for a structure
        repository: (object) The StructureRepository that will be associated
        original_file: (string) Path to the original file (CIF, POSCAR, etc)
        tags: (string or list) Tags that will be associated to that structure
        """
        self.properties = None

        if identifier is None:
            self.structure = structure
            self.identifier = str(_uuid.uuid4())
            if original_file is not None:
                assert (_os.path.isfile(original_file))
            self.original_file = original_file
            self.parents = None
            self.children = None
            if isinstance(tags, str):
                self.tags = [tags]
            elif isinstance(tags, list):
                self.tags = tags
            elif tags is None:
                self.tags = []
            else:
                raise ValueError('The variable tags must be a string or list of strings')

        else:
            assert (original_file is None)
            assert (structure is None)
            assert (tags is None)
            assert (repository is not None)
            self.identifier = identifier
            self.repository = repository
            path = repository.path + '/' + self.identifier
            if not _os.path.isdir(path):
                raise ValueError("Directory not found: " + path)
            if not _os.path.isfile(path + '/metadata.json'):
                raise ValueError("No metadata found in " + path)
            if not _os.path.isfile(path + '/structure.json'):
                raise ValueError("No structure found in " + path)
            self.load()

    def metadatatodict(self):
        ret = {'tags': self.tags,
               'parents': self.parents,
               'children': self.children}
        return ret

    def load(self):
        assert isinstance(self.identifier, str)
        path = self.repository.path + '/' + self.identifier
        rf = open(path + '/metadata.json', 'r')
        self.metadatafromdict(unicode2string(_json.load(rf)))
        rf.close()
        self.structure = load_structure_json(path + '/structure.json')
        if _os.path.isfile(path + '/properties.json'):
            rf = open(path + '/properties.json', 'r')
            self.properties = unicode2string(_json.load(rf))
            rf.close()
        else:
            pass
        if _os.path.isdir(path + '/original') and len(_os.listdir(path + '/original')) > 0:
            self.original_file = path + '/original/' + _os.listdir(path + '/original')[0]

    def save(self):
        path = self.repository.path + '/' + self.identifier
        wf = open(path + '/metadata.json', 'w')
        _json.dump(self.metadatatodict(), wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()
        self.structure.save_json(path + '/structure.json')
        if self.properties is not None:
            wf = open(path + '/properties.json', 'w')
            _json.dump(self.properties, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
        if self.original_file is not None:
            if not _os.path.isdir(path + '/original'):
                _os.mkdir(path + '/original')
            _shutil.copy2(self.original_file, path + '/original')

    def metadatafromdict(self, entrydict):
        self.tags = entrydict['tags']
        self.parents = entrydict['parents']
        self.children = entrydict['children']

    def add_tags(self, tags):
        _add2list(tags, self.tags)

    def add_parents(self, parents):
        _add2list(parents, self.parents)

    def add_children(self, children):
        _add2list(children, self.children)

    def __str__(self):
        ret = 'Structure: \n' + str(self.structure)
        ret += '\nTags: ' + str(self.tags)
        ret += '\nParents: ' + str(self.parents)
        ret += '\nChildren: ' + str(self.children)
        ret += '\nIdentifier: ' + str(self.identifier)
        return ret

    def __eq__(self, other):
        ret = True
        if self.structure != other.structure:
            ret = False
        elif set(self.children) != set(other.children):
            ret = False
        elif set(self.parents) != set(other.parents):
            ret = False
        elif set(self.tags) != set(other.tags):
            ret = False


class ExecutionEntry():
    """
    Defines one execution in the Execution Repository
    """

    def __init__(self, path):
        """
        Creates a new execution repository
        """
        self.path = path


class StructureRepository():
    """
    Defines the location of the executions repository
    and structure repository and methods to add, remove
    and check those repositories
    """

    def __init__(self, path):
        """
        Creates new repositories for calculations and structures

        Args:
        path: (string) Directory path for the structure repository
        """
        self.path = _os.path.abspath(path)

        if _os.path.isfile(self.path + '/repo.json'):
            self.load()
        else:
            self.nentries = 0
            self.tags = {}

            if _os.path.lexists(self.path):
                if not _os.path.isdir(self.path):
                    raise ValueError('Path exists already and it is not a directory')
            else:
                _os.mkdir(self.path)
            self.save()

    def todict(self):
        """
        Serialize the values of the repositories into a dictionary
        """
        repos_dict = {'tags': self.tags,
                      'nentries': self.nentries}

        return repos_dict

    def fromdict(self, repos_dict):

        self.nentries = repos_dict['nentries']
        self.tags = repos_dict['tags']

    def save(self):
        """
        Save an existing repository information
        """
        wf = open(self.path + '/repo.json', 'w')
        _json.dump(self.todict(), wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def load(self):
        """
        Loads an existing repositories from its configuration file
        """
        rf = open(self.path + '/repo.json', 'r')
        self.fromdict(unicode2string(_json.load(rf)))
        rf.close()

    @property
    def get_all_entries(self):
        return [x for x in _os.listdir(self.path) if _os.path.isfile(self.path + '/' + x + '/metadata.json')]

    def merge(self, other):
        """
        Add all the contents from other repositories into the
        calling object

        :param other: StructureRepository
        """
        conflict_entries = []
        for i in other.get_all_enties:
            if i in self.get_all_entries:
                other_structure = StructureEntry(repository=other, identifier=i)
                this_structure = StructureEntry(repository=self, identifier=i)
                if this_structure != other_structure:
                    conflict_entries.append(i)
        if len(conflict_entries) == 0:
            for i in other.get_all_enties:
                if i not in self.get_all_entries:
                    _shutil.copytree(other.path+'/'+i, self.path+'/'+i)
        else:
            print('Conflict entries found, No merge done')
            return conflict_entries

    def add_entry(self, entry):
        """
        Add a new StructureEntry into the repository
        """
        entry.repository = self
        path = self.path + '/' + entry.identifier
        if not _os.path.isdir(path):
            _os.mkdir(path)
        entry.save()
        self.nentries += 1
        if entry.tags is not None:
            for itag in entry.tags:
                if itag in self.tags:
                    if not entry.identifier in self.tags[itag]:
                        self.tags[itag].append(entry.identifier)
                else:
                    self.tags[itag] = [entry.identifier]
        self.save()

    def __str__(self):
        ret = 'Location: ' + self.path
        ret += '\nNumber of entries: ' + str(self.nentries)
        if len(self.tags) > 0:
            for itag in self.tags:
                ret += '\n\t' + itag + ':'
                ret += '\n' + str(self.tags[itag])
        else:
            ret += '\nTags: ' + str(self.tags)
        return ret


class ExecutionRepository():
    """
    Defines the location and properties of the Repository
    where all the executions will be stored
    """

    def __init__(self):
        """
        Creates a Repository for Executions
        """
        pass


def _add2list(orig, dest):
    if isinstance(orig, str):
        if not orig in dest:
            dest.append(orig)
    elif isinstance(orig, list):
        for iorig in dest:
            if not iorig in dest:
                dest.append(iorig)
