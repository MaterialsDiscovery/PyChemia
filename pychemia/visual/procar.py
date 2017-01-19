#!/usr/bin/python
# procar.py
#
# Developers:
#
# -Aldo Romero: hromero@mpi-halle.mpg.de
# -Francisco Munoz: fvmunoz@gmail.com
#


"""
Implemented Classes:

 -UtilsProcar: handy methods, not intented to be used directly by the
  user

 -ProcarParser: reads data from a procar (may be compressed) and store
    them as arrays:
    -kpoint[kpointsCount][3]
    -bands[kpointsCount][bandsCount]
    -spd[kpoint][band][ispin][atom][orbital]

 -ProcarFileFilter: Filter/average a PROCAR, writing a new file with
    the changes. Useful methods
    -FilterOrbitals: write just some orbitals or combination of them
    -FilterAtoms: write just selected atoms (or add them)
    -FilterBands: Writes just selected bands.
    -FilterSpin: Write just selected ispin components (ie: density)

 -ProcarSelect: Select/averages the data yielding a bidimensional
  array ready to plot

"""
# basic modules. Should be present
import logging
import re
import matplotlib.pyplot as plt
import numpy as np


class UtilsProcar:
    """
    This class store handy methods that do not fit any other place

    members:

    -Openfile: Tries to open a File, it has suitable values for PROCARs
     and can handle gzipped files

    -MergeFiles: concatenate two or more PROCAR files taking care of
     metadata and kpoint indexes. Useful for splitted bandstructures
     calculation.

    -FermiOutcar: it greps the Fermi Energy from a given outcar file.

    -RecLatOutcar: it greps the reciprocal lattice from the outcar.

    """

    def __init__(self, loglevel=logging.WARNING):
        self.log = logging.getLogger("UtilsProcar")
        self.log.setLevel(loglevel)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s:"
                                               " %(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)
        self.log.debug("UtilsProcar()")
        self.log.debug("UtilsProcar()...done")
        return

    def open_file(self, filename=None):
        """
        Tries to open a File, it has suitable values for PROCAR and can
        handle gzipped files

        :param filename: (str) Filename to open the PROCAR file
        :return:

        Example:

>>> foo = UtilsProcar.open_file()
Tries to open "PROCAR", then "PROCAR.gz"

>>> foo = UtilsProcar.open_file(filename="../bar")
Tries to open "../bar". If it is a directory, it will try to open
"../bar/PROCAR" and if fails again "../bar/PROCAR.gz"

>>> foo = UtilsProcar.open_file(filename="PROCAR-spd.gz")
Tries to open a gzipped file "PROCAR-spd.gz"

        If unable to open a file, it raises a "IOError" exception.
    """
        import os
        import gzip

        self.log.debug("open_file()")
        self.log.debug("Filename :" + filename)

        if filename is None:
            filename = "PROCAR"
            self.log.debug("Input was None, now is: " + filename)

        # checking if fileName is just a path and needs a "PROCAR to " be
        # appended
        elif os.path.isdir(filename):
            self.log.info("The filename is a directory")
            if filename[-1] != r"/":
                filename += "/"
            filename += "PROCAR"
            self.log.debug("I will try  to open :" + filename)

        # checking that the file exist
        if os.path.isfile(filename):
            self.log.debug("The File does exist")
            # Checking if compressed
            if filename[-2:] == "gz":
                self.log.info("A gzipped file found")
                in_file = gzip.open(filename, "r")
            else:
                self.log.debug("A normal file found")
                in_file = open(filename, "r")
            return in_file

        # otherwise a gzipped version may exist
        elif os.path.isfile(filename + ".gz"):
            self.log.info("File not found, however a .gz version does exist and will"
                          " be used")
            in_file = gzip.open(filename + ".gz")

        else:
            self.log.debug("File not exist, neither a gzipped version")
            raise IOError("File not found")

        self.log.debug("OpenFile()...done")
        return in_file

    def MergeFiles(self, in_files, out_file, gzipOut=False):
        """
        Concatenate two or more PROCAR files. This methods
        takes care of the k-indexes.

        Useful when large number of K points have been calculated in
        different PROCARs.

        Args:
        -inFiles: an iterable with files to be concatenated

        -outFile: a string with the outfile name.

        -gzipOut: whether gzip or not the outout file.

        Warning: spin polarized case is not Ok!

        """
        import gzip

        self.log.debug("MergeFiles()")
        self.log.debug("infiles: " " ,".join(in_files))

        in_files = [self.open_file(x) for x in in_files]
        header = [x.readline() for x in in_files]
        self.log.debug("All the input headers are: \n" + "".join(header))
        metas = [x.readline() for x in in_files]
        self.log.debug("All the input metalines are:\n " + "".join(metas))
        # parsing metalines

        parsedMeta = [map(int, re.findall(r"#[^:]+:([^#]+)", x)) for x in metas]
        kpoints = [x[0] for x in parsedMeta]
        bands = set([x[1] for x in parsedMeta])
        ions = set([x[2] for x in parsedMeta])

        # checking that bands and ions macht (mind: bands & ions are 'sets'):
        if len(bands) != 1 or len(ions) != 1:
            self.log.error("Number of bands/ions  do not match")
            raise RuntimeError("Files are incompatible")

        newKpoints = np.array(kpoints, dtype=int).sum()
        self.log.info("New number of Kpoints: " + str(newKpoints))
        newMeta = metas[0].replace(str(kpoints[0]), str(newKpoints), 1)
        self.log.debug("New meta line:\n" + newMeta)

        if gzipOut:
            self.log.debug("gzipped output")
            out_file = gzip.open(out_file, 'w')
        else:
            self.log.debug("normal output")
            out_file = open(out_file, 'w')
        out_file.write(header[0])
        out_file.write(newMeta)

        # embedded function to change old k-point indexes by the correct
        # ones. The `kreplace.k` syntax is for making the variable 'static'
        def kreplace(matchobj):
            # self.log.debug(print matchobj.group(0))
            kreplace.k += 1
            kreplace.localCounter += 1
            return matchobj.group(0).replace(str(kreplace.localCounter),
                                             str(kreplace.k))

        kreplace.k = 0
        self.log.debug("Going to replace K-points indexes")
        for inFile in in_files:
            lines = inFile.read()
            kreplace.localCounter = 0
            lines = re.sub('(\s+k-point\s*\d+\s*:)', kreplace, lines)
            out_file.write(lines)

        self.log.debug("Closing output file")
        out_file.close()
        self.log.debug("MergeFiles()...done")
        return

    def FermiOutcar(self, filename):
        """Just finds all E-fermi fields in the outcar file and keeps the
        last one (if more than one found).

        :param filename: the file name of the outcar to be readed
        :return:
        """
        self.log.debug("FermiOutcar(): ...")
        self.log.debug("Input filename : " + filename)

        outcar = open(filename, "r").read()
        match = re.findall(r"E-fermi\s*:\s*(-?\d+.\d+)", outcar)[-1]
        self.log.info("Fermi Energy found : " + match)
        self.log.debug("FermiOutcar(): ...Done")
        return float(match)

    def RecLatOutcar(self, filename):
        """Finds and return the reciprocal lattice vectors, if more than
        one set present, it return just the last one.

        Args:
        -filename: the name of the outcar file  to be read

        """
        self.log.debug("RecLatOutcar(): ...")
        self.log.debug("Input filename : " + filename)

        outcar = open(filename, "r").read()
        # just keeping the last component
        recLat = re.findall(r"reciprocal\s*lattice\s*vectors\s*([-.\s\d]*)",
                            outcar)[-1]
        self.log.debug("the match is : " + recLat)
        recLat = recLat.split()
        recLat = np.array(recLat, dtype=float)
        # up to now I have, both direct and rec. lattices (3+3=6 columns)
        recLat.shape = (3, 6)
        recLat = recLat[:, 3:]
        self.log.info("Reciprocal Lattice found :\n" + str(recLat))
        self.log.debug("RecLatOutcar(): ...Done")
        return recLat

    def ProcarRepair(self, infilename, outfilename):
        """It Tries to repair some stupid problems due the stupid fixed
        format of the stupid fortran.

        Up to now it only separes k-points as the following:
        k-point    61 :    0.00000000-0.50000000 0.00000000 ...
        to
        k-point    61 :    0.00000000 -0.50000000 0.00000000 ...

        But as I found new stupid errors they should be fixed here.
        """
        self.log.debug("ProcarRepair(): ...")
        infile = self.open_file(infilename)
        fileStr = infile.read()
        infile.close()

        # Fixing bands issues (when there are more than 999 bands)
        # band *** # energy    6.49554019 # occ.  0.00000000
        fileStr = re.sub(r'(band\s)(\*\*\*)', r'\1 1000', fileStr)

        # Fixing k-point issues
        fileStr = re.sub(r'(\.\d{8})(\d{2}\.)', r'\1 \2', fileStr)
        fileStr = re.sub(r'(\d)-(\d)', r'\1 -\2', fileStr)

        fileStr = re.sub(r'\*+', r' -10.0000 ', fileStr)

        outfile = open(outfilename, 'w')
        outfile.write(fileStr)
        outfile.close()

        self.log.debug("ProcarRepair(): ...Done")
        return


class ProcarParser:
    """
    Parses a PROCAR file and store it in memory. It only deals with
    PROCAR files, that means no Fermi energy (UtilsProcar.FermiOutcar
    can help), and the reciprocal vectors should be supplied (if used,
    see UtilsProcar class).

    Members:

        __init__(self, loglevel): The setup the variables internally, `loglevel`
        sets the verbosity level ie: `loglevel=logging.DEBUG` for debugging. Its
        default is `logging.WARNING`

    Don't use the other methods beggining with underscores "_"

    Example:
    To read a PROCAR or PROCAR.gz file

    >>> foo = ProcarParser()
    >>> foo.readFile()

    #To include the reciprocal vectors, and file name MyFirstPROCAR'

    >>> outcarparser = UtilsProcar()
    >>> recLat = outcarparser.RecLatOutcar(args.outcar)
    >>> foo = ProcarParser()
    >>> foo.readFile("MyFirstPROCAR", recLat=recLat)

    """

    def __init__(self, loglevel=logging.WARNING):
        # array with k-points, they have the following values
        # -None: if not parsed (yet) or parsed with a `permissive` flag on
        # -direct coordinates: if a recLattice was not supplied to the parser
        # -cartesian coords: if a recLattice was supplied to the parser.
        # In the later cases, self.kpoints.shape=(self.kpointsCount, 3)
        self.kpoints = None
        # Number of kpoints, as given by the KPOINTS header (PROCAR file)
        self.kpointsCount = None

        # bands headers present in PROCAR file.
        # self.bands.shape=(self.kpointsCount,self.bandsCount)
        self.bands = None
        # Number of bands. For a spin polarized calculation the number of
        # bands is double (spin ip + spin down). On this array there is no
        # distinction between spin up and down
        self.bandsCount = None

        # Number of ions+1 the +1 is the 'tot' field, ie: the sum over all atoms
        self.ionsCount = None

        self.fileStr = None  # the actual file, stored in memory
        self.spd = None  # the atom/orbital projected data
        self.orbitalName = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz",
                            "dx2", "tot"]
        self.orbitalCount = None  # number of orbitals

        # number of spin components (blocks of data), 1: non-magnetic non
        # polarized, 2: spin polarized collinear, 4: non-collinear
        # spin.
        # NOTE: before calling to `self._readOrbital` the case '4'
        # is marked as '1'
        self.ispin = None
        self.recLattice = None  # reciprocal lattice vectors
        self.utils = UtilsProcar()

        self.log = logging.getLogger("ProcarParser")
        self.log.setLevel(loglevel)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s:"
                                               " %(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)
        # At last, one message to the logger.
        self.log.debug("Procar instanciated")
        return

    def _readKpoints(self, permissive=False):
        """Reads the k-point headers. A typical k-point line is:
        k-point    1 :    0.00000000 0.00000000 0.00000000  weight = 0.00003704\n

        fills self.kpoint[kpointsCount][3]

        The weights are discarded (are they useful?)
        """
        self.log.debug("readKpoints")
        if not self.fileStr:
            self.log.warning("You should invoke `procar.readFile()` instead. Returning")
            return

        # finding all the K-points headers
        self.kpoints = re.findall(r"k-point\s+\d+\s*:\s+([-.\d\s]+)", self.fileStr)
        self.log.debug(str(len(self.kpoints)) + " K-point headers found")
        self.log.debug("The first match found is: " + str(self.kpoints[0]))

        # trying to build an array
        self.kpoints = [x.split() for x in self.kpoints]
        try:
            self.kpoints = np.array(self.kpoints, dtype=float)
        except ValueError:
            self.log.error("Ill-formatted data:")
            print('\n'.join([str(x) for x in self.kpoints]))
            if permissive is True:
                # Discarding the kpoints list, however I need to set
                # self.ispin beforehand.
                if len(self.kpoints) == self.kpointsCount:
                    self.ispin = 1
                elif len(self.kpoints) == 2 * self.kpointsCount:
                    self.ispin = 2
                else:
                    raise ValueError("Kpoints do not match with ispin=1 or 2.")
                self.kpoints = None
                self.log.warning("K-points list is useless, setting it to `None`")
                return
            else:
                raise ValueError("Badly formated Kpoints headers, try `--permissive`")
        # if successful, go on

        # trying to identify an non-polarized or non-collinear case, a
        # polarized case or a defective file

        if len(self.kpoints) != self.kpointsCount:
            # if they do not match, may means two things a spin polarized
            # case or a bad file, lets check
            self.log.debug("Number of kpoints do not match, looking for a "
                           "spin-polarized case")
            # lets start testing if it is spin polarized, if so, there
            # should be 2 identical blocks of kpoints.
            up, down = np.vsplit(self.kpoints, 2)
            if (up == down).all():
                self.log.info("Spin-polarized calculation found")
                self.ispin = 2
                # just keeping one set of kpoints (the other will be
                # discarded)
                self.kpoints = up
            else:
                self.log.error("Number of K-points do not match! check them.")
                raise RuntimeError("Bad Kpoints list.")
        # if ISPIN != 2 setting ISPIN=1 (later for the non-collinear case 1->4)
        # It is unknown until parsing the projected data
        else:
            self.ispin = 1

        # checking again, for compatibility,
        if len(self.kpoints) != self.kpointsCount:
            raise RuntimeError("Kpoints number do not match with metadata (header of PROCAR)")

        self.log.debug(str(self.kpoints))
        self.log.info("The kpoints shape is " + str(self.kpoints.shape))

        if self.recLattice is not None:
            self.log.info("Changing to cartesians coordinates")
            self.kpoints = np.dot(self.kpoints, self.recLattice)
            self.log.debug("New kpoints: \n" + str(self.kpoints))
        return

    def _readBands(self):
        """Reads the bands header. A typical bands is:
        band   1 # energy   -7.11986315 # occ.  1.00000000

        fills self.bands[kpointsCount][bandsCount]

        The occupation numbers are discarded (are they useful?)"""
        self.log.debug("readBands")
        if not self.fileStr:
            self.log.warning("You should invoke `procar.read()` instead. Returning")
            return

        # finding all bands
        self.bands = re.findall(r"band\s*(\d+)\s*#\s*energy\s*([-.\d\s]+)",
                                self.fileStr)
        self.log.debug(str(len(self.bands)) +
                       " bands headers found, bands*Kpoints = " +
                       str(self.bandsCount * self.kpointsCount))
        self.log.debug("The first match found is: " + str(self.bands[0]))

        # checking if the number of bands match

        if len(self.bands) != self.bandsCount * self.kpointsCount * self.ispin:
            self.log.error("Number of bands headers do not match")
            raise RuntimeError("Number of bands don't match")

        # casting to array to manipulate the bands
        self.bands = np.array(self.bands, dtype=float)
        self.log.debug(str(self.bands))

        # Now I will deal with the spin polarized case. The goal is join
        # them like for a non-magnetic case
        if self.ispin == 2:
            # up and down are along the first axis
            up, down = np.vsplit(self.bands, 2)
            self.log.debug("up   , " + str(up.shape))
            self.log.debug("down , " + str(down.shape))

            # reshapping (the 2  means both band index and energy)
            up.shape = (self.kpointsCount, self.bandsCount, 2)
            down.shape = (self.kpointsCount, self.bandsCount, 2)

            # setting the correct number of bands (up+down)
            self.bandsCount *= 2
            self.log.debug("New number of bands : " + str(self.bandsCount))

            # and joining along the second axis (axis=1), ie: bands-like
            self.bands = np.concatenate((up, down), axis=1)

        # otherwise just reshaping is needed
        else:
            self.bands.shape = (self.kpointsCount, self.bandsCount, 2)

        # Making a test if the broadcast is rigth, otherwise just print
        test = [x.max() - x.min() for x in self.bands[:, :, 0].transpose()]
        if np.array(test).any():
            self.log.warning("The indexes of bands do not match. CHECK IT. "
                             "Likely the data was wrongly broadcasted")
            self.log.warning(str(self.bands[:, :, 0]))
        # Now safely removing the band index
        self.bands = self.bands[:, :, 1]
        self.log.info("The bands shape is " + str(self.bands.shape))
        return

    def _readOrbital(self):
        """Reads all the spd-projected data. A typical/expected block is:
        ion      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot
          1  0.079  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.079
          2  0.152  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.152
          3  0.079  0.000  0.001  0.000  0.000  0.000  0.000  0.000  0.000  0.079
          4  0.188  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.188
          5  0.188  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.188
        tot  0.686  0.000  0.002  0.000  0.000  0.000  0.000  0.000  0.000  0.688
        (x2 for spin-polarized -akwardkly formatted-, x4 non-collinear -nicely
        formatted-).

        The data is stored in an array self.spd[kpoint][band][ispin][atom][orbital]

        Undefined behavior in case of phase factors (LORBIT = 12).
        """
        self.log.debug("readOrbital")
        if not self.fileStr:
            self.log.warning("You should invoke `procar.readFile()` instead. Returning")
            return

        # finding all orbital headers
        self.spd = re.findall(r"ion(.+)", self.fileStr)
        self.log.info("the first orbital match reads: " + self.spd[0])
        self.log.debug("And I found " + str(len(self.spd)) + " orbitals headers")

        # testing if the orbital names are known (the standard ones)
        FoundOrbs = self.spd[0].split()
        size = len(FoundOrbs)
        # only the first 'size' orbital
        StdOrbs = self.orbitalName[:size - 1] + self.orbitalName[-1:]
        if FoundOrbs != StdOrbs:
            self.log.warning(str(size) + " orbitals. (Some of) They are unknow (if "
                                         "you did 'filter' them it is OK).")
        self.orbitalCount = size
        self.orbitalNames = self.spd[0].split()
        self.log.debug("Anyway, I will use the following set of orbitals: "
                       + str(self.orbitalNames))

        # Now reading the bulk of data
        self.log.debug("Now searching the values")
        # The case of just one atom is handled differently since the VASP
        # output is a little different
        if self.ionsCount is 1:
            self.spd = re.findall(r"^(\s*1\s+.+)$", self.fileStr, re.MULTILINE)
        else:
            self.spd = re.findall(r"([-.\d\se]+tot.+)\n", self.fileStr)
        # free the memory (could be a lot)
        self.fileStr = None
        self.log.debug("the first entry is \n" + self.spd[0])

        # Now the method will try to find the value of self.ispin,
        # previously it was set to either 1 or 2. If "1", it could be 1 or
        # 4, but previously it was impossible to find the rigth value. If
        # "2" it has to macth with the number of entries of spd data.

        self.log.debug("Number of entries found: " + str(len(self.spd)))
        expected = self.bandsCount * self.kpointsCount
        self.log.debug("The number of entries for a non magnetic calc. is: " +
                       str(expected))
        if expected == len(self.spd):
            self.log.info("Both numbers match, ok, going ahead")
        # catching a non-collinear calc.
        elif expected * 4 == len(self.spd):
            self.log.info("non-collinear calculation found")
            # testing if previous ispin value is ok
            if self.ispin != 1:
                self.log.warning("Incompatible data: self.ispin= " + str(self.ispin) +
                                 ". Now is 4")
            self.ispin = 4
        else:
            self.log.error("The parser or data is wrong!!!")
            self.log.info("bandsCount: " + str(self.bandsCount))
            self.log.info("KpointsCount: " + str(self.kpointsCount))
            raise RuntimeError("Shit happens")

        # checking for consistency
        for line in self.spd:
            if len(line.split()) != self.ionsCount * (self.orbitalCount + 1):
                self.log.error("Expected: " + str(self.ionsCount) + "*" +
                               str(self.orbitalCount + 1) + " = " +
                               str(self.ionsCount * (self.orbitalCount + 1)) +
                               " Fields. Present block: " + str(len(line.split())))
                print(line)
                raise RuntimeError("Flats happens")

        # replacing the "tot" string by a number, to allows a conversion
        # to numpy
        self.spd = [x.replace('tot', '0') for x in self.spd]
        self.spd = [x.split() for x in self.spd]
        self.spd = np.array(self.spd, dtype=float)
        self.log.debug("The spd (old) array shape is:" + str(self.spd.shape))

        # handling collinear polarized case
        if self.ispin == 2:
            self.log.debug("Handling spin-polarized collinear case...")
            # splitting both spin components, now they are along k-points
            # axis (1st axis) but, then should be concatenated along the
            # bands.
            up, down = np.vsplit(self.spd, 2)
            # ispin = 1 for a while, we will made the distinction
            up.shape = (self.kpointsCount, self.bandsCount / 2, 1,
                        self.ionsCount, self.orbitalCount + 1)
            down.shape = (self.kpointsCount, self.bandsCount / 2, 1,
                          self.ionsCount, self.orbitalCount + 1)
            # concatenating bandwise. Density and magntization, their
            # meaning is obvious, and do uses 2 times more memory than
            # required, but I *WANT* to keep it as close as possible to the
            # non-collinear or non-polarized case
            density = np.concatenate((up, down), axis=1)
            magnet = np.concatenate((up, -down), axis=1)
            # concatenated along 'ispin axis'
            self.spd = np.concatenate((density, magnet), axis=2)
            self.log.debug("polarized collinear spd.shape= " + str(self.spd.shape))

        # otherwise, just a reshaping suffices
        else:
            self.spd.shape = (self.kpointsCount, self.bandsCount, self.ispin,
                              self.ionsCount, self.orbitalCount + 1)

        self.log.info("spd array ready. Its shape is:" + str(self.spd.shape))
        return

    def readFile(self, procar=None, permissive=False, recLattice=None):
        """
        Reads and parses the whole PROCAR file. This method is a sort
        of metamethod: it opens the file, reads the meta data and call the
        respective functions for parsing kpoints, bands, and projected
        data.

        The only method of the API it load the file completely.

        :param procar: name of the PROCAR file, can be a gzipped file (the extension is no required).
                        The default covers a wide range of obvious alternatives. The file name, if `None`
                        or a directory, a suitable set of defaults will be used. Default=None

        :param permissive: Set to `True` if the PROCAR file has problems reading
                            the Kpoints (stupid Fortran), but in that case the
                            Kpoints mesh will be discarded. Future updates could
                            allow it to handle other formating/corruption issues.
                            turn on (or off) some features to deal with badly
                            written PROCAR files (stupid fortran), up to now just ignores the
                            kpoints coordinates, which -as side effect- prevent he rigth
                            space between kpoints. Default=False (off)

        :param recLattice: Reciprical Vectors, you want to provide them since not
                            all the paths on the BZ are the same.  a 3x3 array containing the reciprocal vectors, to
                            change the Kpoints from rec. coordinates to cartesians. Rarely
                            given by hand, see `UtilsProcar.RecLatProcar`. If given, the
                            kpoints will be converted from direct coordinates to cartesian
                            ones. Default=None
        :return:

        """
        self.log.debug("readFile...")

        self.recLattice = recLattice

        self.log.debug("Opening file: '" + str(procar) + "'")
        f = self.utils.open_file(procar)
        # Line 1: PROCAR lm decomposed
        f.readline()  # throwaway
        # Line 2: # of k-points:  816   # of bands:  52   # of ions:   8
        metaLine = f.readline()  # metadata
        self.log.debug("The metadata line is: " + metaLine)
        re.findall(r"#[^:]+:([^#]+)", metaLine)
        self.kpointsCount, self.bandsCount, self.ionsCount = \
            map(int, re.findall(r"#[^:]+:([^#]+)", metaLine))
        self.log.info("kpointsCount = " + str(self.kpointsCount))
        self.log.info("bandsCount = " + str(self.bandsCount))
        self.log.info("ionsCount = " + str(self.ionsCount))
        if self.ionsCount is 1:
            self.log.warning("Special case: only one atom found. The program may not work as expected")
        else:
            self.log.debug("An extra ion representing the  total value will be added")
            self.ionsCount += 1

        # reading all the rest of the file to be parsed below
        self.fileStr = f.read()
        self._readKpoints(permissive)
        self._readBands()
        self._readOrbital()
        self.log.debug("readfile...done")
        return


class ProcarFileFilter:
    """Process a PROCAR file fields line-wise, specially useful for HUGE
    files. This could be thought as pre-processing, writting a new
    PROCAR-like file but reduced in some way.

    A PROCAR File is basically an multi-dimmensional arrays of data, the
    largest being:
    spd_data[#kpoints][#band][#ispin][#atom][#orbital]

    while the number of Kpoints d'ont seems a target for reduction
    (omission or averaging), the other fields can be reduced, for
    instance: grouping the atoms by species or as "surtrate" and
    "adsorbate", or just keeping the bands close to the Fermi energy, or
    discarding the d-orbitals in a s-p system. You got the idea, rigth?

    Example:

    -To group the "s", "p" y "d" orbitals from the file PROCAR and write
     them in PROCAR-spd:

     >>> a = ProcarFileFilter("PROCAR", "PROCAR-new")
     >>> a.FilterOrbitals([[0], [1, 2, 3], [4, 5, 6, 7, 8]], ['s', 'p', 'd'])

         The PROCAR-new will have just 3+1 columns (orbitals are colum-wise
         , take a look to the file). If you omit the ['s', 'p', 'd'] list,
         the new orbitals will have a generic meaningless name (o1, o2, o3)

    -To group the atoms 1,2,5,6 and 3,4,7,8 from PROCAR and write them
     in PROCAR-new (note the 0-based indexes):

     >>> a = ProcarFileFilter("PROCAR", "PROCAR-new")
     >>> a.FilterAtoms([[0, 1, 4, 5], [2, 3, 6, 7]])

     -To select just the total density (ie: ignoring the spin-resolved stuff,
      if any)from PROCAR and write it in PROCAR-new:

     >>> a = ProcarFileFilter("PROCAR", "PROCAR-new")
     >>> a.FilterSpin([0])

    """

    def __init__(self, infile=None, outfile=None, loglevel=logging.WARNING):
        """Initialize the class.

        Params: `infile=None`, input fileName
        """
        self.infile = infile
        self.outfile = outfile

        # We want a logging to tell us what is happening
        self.log = logging.getLogger("ProcarFileFilter")
        self.log.setLevel(loglevel)
        # This is a handler for logging, by now just keep it
        # untouched. Dont really matters its usage
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s:"
                                               " %(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)
        # At last, one message to the logger.
        self.log.debug("ProcarFileFilter instanciated")
        return

    def setInFile(self, infile):
        """Sets a input file `infile`, it can contains the path to the file"""
        self.infile = infile
        self.log.info("Input File: " + infile)
        return

    def setOutFile(self, outfile):
        """Sets a output file `outfile`, it can contains the path to the file"""
        self.outfile = outfile
        self.log.info("Out File: " + outfile)
        return

    def FilterOrbitals(self, orbitals, orbitalsNames):
        """
        Reads the file already set by SetInFile() and writes a new
        file already set by SetOutFile(). The new file only has the
        selected/grouped orbitals.

        :param orbitals: nested iterable with the orbitals indexes to be
            considered. For example: [[0],[2]] means select the first
            orbital ("s") and the second one ("pz").
            [[0],[1,2,3],[4,5,6,7,8]] is ["s", "p", "d"].
        :param orbitalsNames: The name to be put in each new orbital field (of a
            orbital line). For example ["s","p","d"] is a good
            orbitalsName for the orbitals=[[0],[1,2,3],[4,5,6,7,8]].
            However, ["foo", "bar", "baz"] is equally valid.

        :return:

        Note:
            - The atom index is not counted as the first field.
            - The last column ('tot') is so important that it is always
                included. Do not needs to be called

        """
        # setting iostuff, this method -and class- should not made any
        # checking about IO, that is the job of the caller
        self.log.info("In File: " + self.infile)
        self.log.info("Out File: " + self.outfile)
        # open the files
        fout = open(self.outfile, 'w')
        fopener = UtilsProcar()
        fin = fopener.open_file(self.infile)
        for line in fin:
            if re.match(r"\s*ion\s*", line):
                # self.log.debug("orbital line found: " + line)
                line = " ".join(['ion'] + orbitalsNames + ['tot']) + "\n"

            elif re.match(r"\s*\d+\s*", line) or re.match(r"\s*tot\s*", line):
                # self.log.debug("data line found: " + line)
                line = line.split()
                # all floats to an array
                data = np.array(line[1:], dtype=float)
                # setting a new line, keeping just the first value
                line = line[:1]
                for orbset in orbitals:
                    line.append(data[orbset].sum())
                # the last value ("tot") always  should be written
                line.append(data[-1])
                # converting to str
                line = [str(x) for x in line]
                line = " ".join(line) + "\n"
            fout.write(line)

        return

    def FilterAtoms(self, atomsGroups):
        """
        Reads the file already set by SetInFile() and writes a new
        file already set by SetOutFile(). The new file only has the
        selected/grouped atoms.

        Args:

        -atomsGroups: nested iterable with the atoms indexes (0-based) to
          be considered. For example: [[0],[2]] means select the first and
          the second atoms. While [[1,2,3],[4,5,6,7,8]] means select the
          contribution of atoms 1+2+3 and 4+5+6+7+8

        Note:
          -The atom index is c-based (or python) beginning with 0
          -The output has a dummy atom index, without any intrisic meaning

        """
        # setting iostuff, this method -and class- should not made any
        # checking about IO, that is the job of the caller
        self.log.info("In File: " + self.infile)
        self.log.info("Out File: " + self.outfile)
        # open the files
        fout = open(self.outfile, 'w')
        fopener = UtilsProcar()
        with fopener.open_file(self.infile) as fin:
            # I need to change the numbers of ions, it will needs the second
            # line. The first one is not needed
            fout.write(fin.readline())
            line = fin.readline()
            line = line.split()
            # the very last value needs to be changed
            line[-1] = str(len(atomsGroups))
            line = ' '.join(line)
            fout.write(line + '\n')

            # now parsing the rest of the file
            data = []
            for line in fin:
                # if line has data just capture it
                if re.match(r"\s*\d+\s*", line):
                    # self.log.debug("atoms line found: " + line)
                    data.append(line)
                # if `line` is a end of th block (begins with 'tot'), do the
                # work. And clean up data then
                elif re.match(r"\s*tot\s*", line):
                    # self.log.debug("tot line found: " + line)
                    # making an array
                    data = [x.split() for x in data]
                    data = np.array(data, dtype=float)
                    # iterating on the atoms groups
                    for index in range(len(atomsGroups)):
                        atoms = atomsGroups[index]
                        # summing colum-wise
                        atomLine = data[atoms].sum(axis=0)
                        atomLine = [str(x) for x in atomLine]
                        # the atom index should not be averaged (anyway now is
                        # meaningless)
                        atomLine[0] = str(index + 1)
                        atomLine = ' '.join(atomLine)
                        fout.write(atomLine + '\n')

                    # clean the buffer
                    data = []
                    # and write the `tot` line
                    fout.write(line)
                # otherwise just write this line
                else:
                    fout.write(line)

        return

    def FilterBands(self, Min, Max):
        """
        Reads the file already set by SetInFile() and writes a new
        file already set by SetOutFile(). The new file only has the
        selected bands.

        Args:

        -Min, Max:
          the minimum/maximum band  index to be considered, the indexes
          are the same used by vasp (ie written in the file).


        Note: -Since bands are somewhat disordered in vasp you may like to
          consider a large region and made some trial and error

        """
        # setting iostuff, this method -and class- should not made any
        # checking about IO, that is the job of the caller
        self.log.info("In File: " + self.infile)
        self.log.info("Out File: " + self.outfile)
        # open the files
        fout = open(self.outfile, 'w')
        fopener = UtilsProcar()
        fin = fopener.open_file(self.infile)

        # I need to change the numbers of kpoints, it will needs the second
        # line. The first one is not needed
        fout.write(fin.readline())
        line = fin.readline()
        # the third value needs to be changed, however better print it
        self.log.debug("The line contaning bands number is " + line)
        line = line.split()
        self.log.debug("The number of bands is: " + line[7])
        line[7] = str(Max - Min + 1)
        line = ' '.join(line)
        fout.write(line + '\n')

        # now parsing the rest of the file
        write = True
        for line in fin:
            if re.match(r"\s*band\s*", line):
                # self.log.debug("bands line found: " + line)
                band = int(re.match(r"\s*band\s*(\d+)", line).group(1))
                if band < Min or band > Max:
                    write = False
                else:
                    write = True
            if re.match(r"\s*k-point\s*", line):
                write = True
            if write:
                fout.write(line)
        return

    def FilterSpin(self, components):
        """Reads the file already set by SetInFile() and writes a new
        file already set by SetOutFile(). The new file only has the
        selected part of the density (sigma_i).

        Args:

        -components: The spin component block, for instante [0] menas just
          the density, while [1,2] would be the the sigma_x and sigma_y
          for a non-collinear calculation.

        TODO: spin-polarized collinear case is not included at all, not
        even a warning message!

        """
        # setting iostuff, this method -and class- should not made any
        # checking about IO, that is the job of the caller
        self.log.info("In File: " + self.infile)
        self.log.info("Out File: " + self.outfile)
        # open the files
        fout = open(self.outfile, 'w')
        fopener = UtilsProcar()
        with fopener.open_file(self.infile) as fin:
            counter = 0
            for line in fin:
                # if any data found
                if re.match(r"\s*\d", line):
                    # check if should be written
                    if counter in components:
                        fout.write(line)
                elif re.match(r"\s*tot", line):
                    if counter in components:
                        fout.write(line)
                    # the next block will belong to other component
                    counter += 1
                elif re.match(r"\s*ion", line):
                    fout.write(line)
                    counter = 0
                else:
                    fout.write(line)
        return


class ProcarSelect:
    """
    Reduces the dimensionality of the data making it uselful to
    plot bands.

    The main data to manipulate is the projected electronic structure.
    Its shape original is:

    spd[kpoint][band][ispin][atom][orbital].

    The selection of components should be done in order, says, first
    "ispin", then "atom", and at last "orbital".

    Note: once any selection has been performed, the data itself
    changes. Say, if you want compare atom [0] and [1,2], you need two
    instances of this class.


    Example to compare the bandstructure of two set of atoms
    >>>

    """

    def __init__(self, ProcarData=None, deepCopy=True, loglevel=logging.WARNING):

        self.spd = None
        self.bands = None
        self.kpoints = None

        # We want a logging to tell us what is happening
        self.log = logging.getLogger("ProcarSelect")
        self.log.setLevel(loglevel)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s:"
                                               " %(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)
        # At last, one message to the logger.
        self.log.debug("ProcarSelect: instanciated")

        if ProcarData is not None:
            self.setData(ProcarData, deepCopy)
        return

    def setData(self, ProcarData, deepCopy=True):
        """
        The data from ProcarData is deepCopy-ed by default (ie: their
        elements are not modified by this class.

        Args:

        -ProcarData: is a ProcarParser instance (or anything with similar
         functionality, duck typing)

        -deepCopy=True: If false a shallow copy will be made (saves memory).
        """
        self.log.debug("setData: ...")
        if deepCopy is True:
            self.spd = ProcarData.spd.copy()
            self.bands = ProcarData.bands.copy()
            self.kpoints = ProcarData.kpoints.copy()
        else:
            self.spd = ProcarData.spd
            self.bands = ProcarData.bands
            self.kpoints = ProcarData.kpoints
        self.log.debug("setData: ... Done")
        return

    def selectIspin(self, value):
        """
        value is a list with the values of Ispin to select.

        Examples:
        >>> foo = ProcarParser()
        >>> foo.readFile("PROCAR")
        >>> bar = ProcarSelect(foo)
        >>> bar.selectIspin([0])  # just the density

        """
        # all kpoint, all bands, VALUE spin, all the rest
        self.log.debug("selectIspin: ...")
        self.log.debug("old spd shape =" + str(self.spd.shape))
        # first, testing if the domensionaluty is rigth:
        dimen = len(self.spd.shape)
        if dimen != 5:
            self.log.error("The array is " + str(dimen) + " dimensional, expecting a"
                                                          " 5 dimensional array.")
            self.log.error("You should call selectIspin->selecAtom->selectOrbitals, "
                           "in this order.")
            raise RuntimeError('Wrong dimensionality of the array')
        self.log.debug("ispin value = " + str(value))
        self.spd = self.spd[:, :, value]
        self.spd = self.spd.sum(axis=2)
        self.log.info("new spd shape =" + str(self.spd.shape))
        self.log.debug("selectIspin: ...Done")
        return

    def selectAtoms(self, value, fortran=False):
        """
        value is a list with the values of Atoms to select. The optional
        `fortran` argument indicates whether a c-like 0-based indexing
        (`=False`, default) or a fortran-like 1-based (`=True`) is
        provided in `value`.

        Examples:
        >>> foo = ProcarParser()
        >>> foo.readFile("PROCAR")
        >>> bar = ProcarSelect(foo)
        >>> bar.selectIspin([...])
        >>> bar.selectAtoms([0, 1, 2])  # atom0+atom1+atom2

        Note: this method should be called after select.Ispin
        """
        self.log.debug("selectAtoms: ...")

        # taking care about stupid fortran indexing
        if fortran is True:
            value = [x - 1 for x in value]

        # all kpoint, all bands, VALUE atoms, all the rest
        self.log.debug("old shape =" + str(self.spd.shape))

        # testing if the dimensionaluty is rigth:
        dimen = len(self.spd.shape)
        if dimen != 4:
            self.log.error("The array is " + str(dimen) + " dimensional, expecting a"
                                                          " 4 dimensional array.")
            self.log.error("You should call selectIspin->selecAtom->selectOrbitals, "
                           "in this order.")
            raise RuntimeError('Wrong dimensionality of the array')
        self.spd = self.spd[:, :, value]
        self.spd = self.spd.sum(axis=2)
        self.log.info("new shape =" + str(self.spd.shape))
        self.log.debug("selectAtoms: ...Done")
        return

    def selectOrbital(self, value):
        """
        value is a list with the values of orbital to select.

        Examples:
        >>> foo = ProcarParser()
        >>> foo.readFile("PROCAR")
        >>> bar = ProcarSelect(foo)
        >>> bar.selectIspin([...])
        >>> bar.selectAtoms([...])
        >>> bar.selectOrbital([-1])  # the last (`tot`) field

        to select "p" orbitals just change the argument in the last line
        to [2,3,4] or as needed

        Note: this method should be called after `select.Ispin` and
        `select.Atoms`

        """
        self.log.debug("selectOrbital: ...")
        self.log.debug("Changing the orbital `values` to have a 0-based indexes")
        # Mind: the first orbital field is the atoms number, which is not
        # an orbital, therefore the orbital index is an affective 1-based
        # therefore all `value` indexes += 1 (well, negative values do not
        # change )
        for i in range(len(value)):
            if value[i] >= 0:
                value[i] += 1

        self.log.debug("New values (indexes to select) :" + str(value))

        # all kpoint, all bands, VALUE orbitals, nothing else?
        self.spd = self.spd[:, :, value]
        self.log.debug("old shape =" + str(self.spd.shape))

        # testing if the dimensionaluty is rigth:
        dimen = len(self.spd.shape)
        if dimen != 3:
            self.log.error("The array is " + str(dimen) + " dimensional, expecting a"
                                                          " 3 dimensional array.")
            self.log.error("You should call selectIspin->selecAtom->selectOrbitals, "
                           "in this order.")
            raise RuntimeError('Wrong dimensionality of the array')

        self.spd = self.spd.sum(axis=2)
        self.log.info("new shape =" + str(self.spd.shape))
        self.log.debug("selectOrbital: ...Done")
        return


class ProcarPlot:
    def __init__(self, bands, spd, kpoints=None):
        self.bands = bands.transpose()
        self.spd = spd.transpose()
        self.kpoints = kpoints
        return

    def plotBands(self, size=None, marker='o', ticks=None):
        if size is not None:
            size /= 2
        if self.kpoints is not None:
            xaxis = [0]
            for i in range(1, len(self.kpoints)):
                d = self.kpoints[i - 1] - self.kpoints[i]
                d = np.sqrt(np.dot(d, d))
                xaxis.append(d + xaxis[-1])
            xaxis = np.array(xaxis)
        else:
            xaxis = np.arange(len(self.bands))
        print("self.kpoints: ", self.kpoints.shape)
        print("xaxis.shape : ", xaxis.shape)
        print("bands.shape : ", self.bands.shape)
        plot = plt.plot(xaxis, self.bands.transpose(), 'r-', marker=marker,
                        markersize=size)
        plt.xlim(xaxis.min(), xaxis.max())

        # handling ticks
        if ticks:
            ticks, ticksNames = zip(*ticks)
            ticks = [xaxis[x] for x in ticks]
            plt.xticks(ticks, ticksNames)

        return plot

    def scatterPlot(self, size=50, mask=None, cmap='hot_r', vmax=None, vmin=None,
                    marker='o', ticks=None):
        bsize, ksize = self.bands.shape
        print(bsize, ksize)

        if self.kpoints is not None:
            xaxis = [0]
            for i in range(1, len(self.kpoints)):
                d = self.kpoints[i - 1] - self.kpoints[i]
                d = np.sqrt(np.dot(d, d))
                xaxis.append(d + xaxis[-1])
            xaxis = np.array(xaxis)
        else:
            xaxis = np.arange(ksize)

        xaxis.shape = (1, ksize)
        xaxis = xaxis.repeat(bsize, axis=0)
        if mask is not None:
            mbands = np.ma.masked_array(self.bands, np.abs(self.spd) < mask)
        else:
            mbands = self.bands
        fig, ax = plt.subplots()
        pc = ax.scatter(xaxis, mbands, c=self.spd, s=size, linewidths=0,
                        cmap=cmap, vmax=vmax, vmin=vmin, marker=marker,
                        edgecolors='none')
        fig.colorbar(pc)
        ax.set_xlim(xaxis.min(), xaxis.max())

        # handling ticks
        if ticks:
            ticks, ticksNames = zip(*ticks)
            ticks = [xaxis[0, x] for x in ticks]
            ax.set_xticks(ticks, ticksNames)

        return fig, ax

    def parametricPlot(self, cmap='hot_r', vmin=None, vmax=None, mask=None,
                       ticks=None):
        from matplotlib.collections import LineCollection
        import matplotlib

        fig = plt.figure()
        gca = fig.gca()
        bsize, ksize = self.bands.shape

        # print self.bands
        if mask is not None:
            mbands = np.ma.masked_array(self.bands, np.abs(self.spd) < mask)
        else:
            # Faking a mask, all elemtnet are included
            mbands = np.ma.masked_array(self.bands, False)
        # print mbands

        if vmin is None:
            vmin = self.spd.min()
        if vmax is None:
            vmax = self.spd.max()
        print("normalizing to: ", (vmin, vmax))
        norm = matplotlib.colors.Normalize(vmin, vmax)

        if self.kpoints is not None:
            xaxis = [0]
            for i in range(1, len(self.kpoints)):
                d = self.kpoints[i - 1] - self.kpoints[i]
                d = np.sqrt(np.dot(d, d))
                xaxis.append(d + xaxis[-1])
            xaxis = np.array(xaxis)
        else:
            xaxis = np.arange(ksize)

        for y, z in zip(mbands, self.spd):
            # print xaxis.shape, y.shape, z.shape
            points = np.array([xaxis, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=plt.get_cmap(cmap), norm=norm,
                                alpha=0.8)
            lc.set_array(z)
            lc.set_linewidth(2)
            gca.add_collection(lc)
        plt.colorbar(lc)
        plt.xlim(xaxis.min(), xaxis.max())
        plt.ylim(mbands.min(), mbands.max())

        # handling ticks
        if ticks:
            ticks, ticksNames = zip(*ticks)
            ticks = [xaxis[x] for x in ticks]
            plt.xticks(ticks, ticksNames)

        return fig

    def atomicPlot(self, cmap='hot_r', vmin=None, vmax=None):
        """
        Just a handler to parametricPlot. Useful to plot energy levels.

        It adds a fake k-point. Shouldn't be invoked with more than one
        k-point
        """

        print("Atomic plot: bands.shape  :", self.bands.shape)
        print("Atomic plot: spd.shape    :", self.spd.shape)
        print("Atomic plot: kpoints.shape:", self.kpoints.shape)

        self.bands = np.hstack((self.bands, self.bands))
        self.spd = np.hstack((self.spd, self.spd))
        self.kpoints = np.vstack((self.kpoints, self.kpoints))
        self.kpoints[0][-1] += 1
        print("Atomic plot: bands.shape  :", self.bands.shape)
        print("Atomic plot: spd.shape    :", self.spd.shape)
        print("Atomic plot: kpoints.shape:", self.kpoints.shape)

        print(self.kpoints)

        fig = self.parametricPlot(cmap, vmin, vmax)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())

        # labels on each band
        for i in range(len(self.bands[:, 0])):
            # print i, self.bands[i]
            plt.text(0, self.bands[i, 0], str(i + 1), fontsize=15)

        return fig


class FermiSurface:
    def __init__(self, kpoints, bands, spd, loglevel=logging.WARNING):
        """FermiSurface: Class to build and to plot a 2D fermi surface.  It
        finds the relevant bands (crossig the Fermi level) and interpolate
        them

        args:
          kpoints: Numpy array with kpoints Nx3.
          bands: the bands with Fermi energy already substracted!, numpy array,
                 Nkpoints x Nbands.
          spd: character (atomic, orbital) of each bands at each Kpoint, numpy
               array Nkpoints x Nbands.
          loglevel(=logging.WARNING): the verbosity level.
        """
        # Since some time ago Kpoints are in cartesian coords (ready to use)
        self.kpoints = kpoints
        self.bands = bands
        self.spd = spd
        self.useful = None  # List of useful bands (filled in findEnergy)
        self.energy = None

        self.log = logging.getLogger("FermiSurface")
        self.log.setLevel(loglevel)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s: "
                                               "%(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)

        self.log.debug("FermiSurface.init: ...")
        self.log.info("Kpoints.shape : " + str(self.kpoints.shape))
        self.log.info("bands.shape   : " + str(self.bands.shape))
        self.log.info("spd.shape     : " + str(self.spd.shape))
        self.log.debug("FermiSurface.init: ...Done")
        return

    def FindEnergy(self, energy):
        self.log.debug("FindEnergy: ...")
        self.energy = energy
        self.log.info("Energy   : " + str(energy))
        bands = self.bands.transpose()
        # searching for bands crossing the desired energy
        self.useful = np.where(np.logical_and(bands.min(axis=1) < energy,
                                              bands.max(axis=1) > energy))
        self.log.info("set of useful bands    : " + str(self.useful))
        bands = bands[self.useful]
        self.log.debug("new bands.shape : " + str(bands.shape))
        if len(bands) == 0:
            self.log.error("No bands found in that range. Check your data. Returning")
            raise RuntimeError("No bands to plot")
        self.log.debug("FindEnergy: ...Done")
        return

    def Plot(self, interpolation=200, mask=None):
        """Only 2D layer geometry along z"""
        self.log.debug("Plot: ...")
        from scipy.interpolate import griddata

        if self.useful is None:
            raise RuntimeError("self.FindEnergy() must be called before Plotting")

        # selecting components of K-points
        x, y = self.kpoints[:, 0], self.kpoints[:, 1]
        self.log.debug("k_x[:10], k_y[:10] values" + str([x[:10], y[:10]]))

        bands = self.bands.transpose()[self.useful]

        # and new, interpolated component
        xmax, xmin = x.max(), x.min()
        ymax, ymin = y.max(), y.min()
        self.log.debug("xlim = " + str([xmin, xmax]) + "  ylim = " + str([ymin, ymax]))
        xnew, ynew = np.mgrid[xmin:xmax:interpolation * 1j, ymin:ymax:interpolation * 1j]

        # interpolation
        bnew = []
        for band in bands:
            self.log.debug("Interpolating ...")
            bnew.append(griddata((x, y), band, (xnew, ynew), method='cubic'))

        plots = [plt.contour(xnew, ynew, z, [self.energy], linewidths=0.5, colors='k', linestyles='solid') for z in
                 bnew]
        plt.axis("equal")
        self.log.debug("Plot: ...Done")
        return plots

    def st(self, sx, sy, sz, spin=None, noarrow=False, interpolation=300):
        """Only 2D layer geometry along z. It is like a enhanced version
        of 'plot' method.

        sx, sy, sz are spin projected Nkpoints x Nbands numpy arrays. They
        also are (already) projected by orbital and atom (from other
        class)

        """
        self.log.debug("st: ...")
        from scipy.interpolate import griddata

        if self.useful is None:
            raise RuntimeError("self.FindEnergy() must be called before Plotting")

        # selecting components of K-points
        x, y = self.kpoints[:, 0], self.kpoints[:, 1]

        bands = self.bands.transpose()[self.useful]

        sx = sx.transpose()[self.useful]
        sy = sy.transpose()[self.useful]
        sz = sz.transpose()[self.useful]

        # and new, interpolated component
        xmax, xmin = x.max(), x.min()
        ymax, ymin = y.max(), y.min()
        self.log.debug("xlim = " + str([xmin, xmax]) + "  ylim = " + str([ymin, ymax]))
        xnew, ynew = np.mgrid[xmin:xmax:interpolation * 1j, ymin:ymax:interpolation * 1j]

        # interpolation
        bnew = []
        for band in bands:
            self.log.debug("Interpolating ...")
            bnew.append(griddata((x, y), band, (xnew, ynew), method='cubic'))

        linewidths = 0.7
        if noarrow:
            linewidths = 0.2
        cont = [plt.contour(xnew, ynew, z, [self.energy],
                            linewidths=linewidths, colors='k',
                            linestyles='solid') for z in bnew]
        plt.axis("equal")

        for (contour, spinX, spinY, spinZ) in zip(cont, sx, sy, sz):
            # The previous interp. yields the level curves, nothing more is
            # useful from there
            paths = contour.collections[0].get_paths()
            verts = [xx.vertices for xx in paths]
            points = np.concatenate(verts)
            self.log.debug("Fermi surf. points.shape: " + str(points.shape))

            newSx = griddata((x, y), spinX, (points[:, 0], points[:, 1]))
            newSy = griddata((x, y), spinY, (points[:, 0], points[:, 1]))
            newSz = griddata((x, y), spinZ, (points[:, 0], points[:, 1]))

            self.log.info("newSx.shape: " + str(newSx.shape))

            import matplotlib.colors as colors

            if noarrow is False:
                plt.quiver(points[::6, 0], points[::6, 1], newSx[::6],
                           newSy[::6], newSz[::6], scale_units='xy',
                           angles='xy', norm=colors.Normalize(-0.5, 0.5))
            else:
                # a dictionary to select the right spin component
                spinDict = {1: newSx[::6], 2: newSy[::6], 3: newSz[::6]}
                plt.scatter(points[::6, 0], points[::6, 1], c=spinDict[spin],
                            s=50, edgecolor='none', alpha=1.0, marker=".",
                            cmap='seismic', norm=colors.normalize(-0.5, 0.5))
        plt.colorbar()
        plt.axis("equal")
        font = {'size': 16}
        plt.rc('font', **font)

        self.log.debug("st: ...Done")
        return


def _q_mult(q1, q2):
    """
    Multiplication of quaternions, it doesn't fit in any other place
    """
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return np.array((w, x, y, z))


class ProcarSymmetry:
    def __init__(self, kpoints, bands, character=None, sx=None, sy=None,
                 sz=None, loglevel=logging.WARNING):
        self.log = logging.getLogger("ProcarSymmetry")
        self.log.setLevel(loglevel)
        self.ch = logging.StreamHandler()
        self.ch.setFormatter(logging.Formatter("%(name)s::%(levelname)s: "
                                               "%(message)s"))
        self.ch.setLevel(logging.DEBUG)
        self.log.addHandler(self.ch)
        self.log.debug("ProcarSymmetry.__init__: ...")

        self.kpoints = kpoints
        self.bands = bands
        # optional arguments when not given will False, but they can still
        # be treated like arrays
        self.character = np.array([])
        if character is not None:
            self.character = character
        self.sx = np.array([])
        if sx is not None:
            self.sx = sx
        self.sy = np.array([])
        if sy is not None:
            self.sy = sy
        self.sz = np.array([])
        if sz is not None:
            self.sz = sz

        self.log.info("Kpoints : " + str(self.kpoints.shape))
        self.log.info("bands   : " + str(self.bands.shape))
        self.log.info("character  : " + str(self.character.shape))
        self.log.info("sx      : " + str(self.sx.shape))
        self.log.info("sy      : " + str(self.sy.shape))
        self.log.info("sz      : " + str(self.sz.shape))
        self.log.debug("ProcarSymmetry.__init__: ...Done")

        return

    def GeneralRotation(self, angle, rotAxis=None, store=True):
        """
        Apply a rotation defined by an angle and an axis.

        Returning value: (Kpoints, sx,sy,sz), the rotated Kpoints and spin
                         vectors (if not the case, they will be empty
                         arrays).

        Arguments
        angle: the rotation angle, must be in degrees!

        rotAxis : a fixed Axis when applying the symmetry, usually it is
        from Gamma to another point). It doesn't need to be normalized.
        The RotAxis can be:
        [x,y,z] : a cartesian vector in k-space.
        'x': [1,0,0], a rotation in the yz plane.
        'y': [0,1,0], a rotation in the zx plane.
        'z': [0,0,1], a rotation in the xy plane

        """
        if rotAxis is None:
            rotAxis = [0, 0, 1]
        if rotAxis == 'x' or rotAxis == 'X':
            rotAxis = [1, 0, 0]
        if rotAxis == 'y' or rotAxis == 'Y':
            rotAxis = [0, 1, 0]
        if rotAxis == 'z' or rotAxis == 'Z':
            rotAxis = [0, 0, 1]
        rotAxis = np.array(rotAxis, dtype=float)
        self.log.debug("rotAxis : " + str(rotAxis))
        rotAxis /= np.linalg.norm(rotAxis)
        self.log.debug("rotAxis Normalized : " + str(rotAxis))
        self.log.debug("Angle : " + str(angle))
        angle = angle * np.pi / 180
        # defining a quaternion for rotatoin
        angle /= 2
        rotAxis *= np.sin(angle)
        qRot = np.array((np.cos(angle), rotAxis[0], rotAxis[1], rotAxis[2]))
        qRotI = np.array((np.cos(angle), -rotAxis[0], -rotAxis[1], -rotAxis[2]))
        self.log.debug("Rot. quaternion : " + str(qRot))
        self.log.debug("Rot. quaternion conjugate : " + str(qRotI))
        # converting self.kpoints into quaternions
        w = np.zeros((len(self.kpoints), 1))
        qvectors = np.column_stack((w, self.kpoints)).transpose()
        self.log.debug("Kpoints-> quaternions (transposed):\n" + str(qvectors.transpose()))
        qvectors = _q_mult(qRot, qvectors)
        qvectors = _q_mult(qvectors, qRotI).transpose()
        kpoints = qvectors[:, 1:]
        self.log.debug("Rotated kpoints :\n" + str(qvectors))

        # rotating the spin vector (if exist)
        sxShape, syShape, szShape = self.sx.shape, self.sy.shape, self.sz.shape
        self.log.debug("Spin vector Shapes : " + str((sxShape, syShape, szShape)))
        # The first entry has to be an array of 0s, w could do the work,
        # but if len(self.sx)==0 qvectors will have a non-defined length
        qvectors = (0 * self.sx.flatten(), self.sx.flatten(),
                    self.sy.flatten(), self.sz.flatten())
        self.log.debug("Spin vector quaternions: \n" + str(qvectors))
        qvectors = _q_mult(qRot, qvectors)
        qvectors = _q_mult(qvectors, qRotI)
        self.log.debug("Spin quaternions after rotation:\n" + str(qvectors))
        sx, sy, sz = qvectors[1], qvectors[2], qvectors[3]
        sx.shape, sy.shape, sz.shape = sxShape, syShape, szShape

        if store is True:
            self.kpoints, self.sx, self.sy, self.sz = kpoints, sx, sy, sz
        self.log.debug("GeneralRotation: ...Done")
        return kpoints, sx, sy, sz

    def RotSymmetryZ(self, order):
        """Applies the given rotational crystal symmetry to the current
        system. ie: to unfold the irreductible BZ to the full BZ.

        Only rotations along z-axis are performed, you can use
        self.GeneralRotation first.

        The user is responsible of provide a useful input. The method
        doesn't check the physics.

        """
        self.log.debug("RotSymmetryZ:...")
        rotations = [self.GeneralRotation(360 * i / order, store=False) for i
                     in range(order)]
        rotations = zip(*rotations)
        self.log.debug("self.kpoints.shape (before concat.): " +
                       str(self.kpoints.shape))
        self.kpoints = np.concatenate(rotations[0], axis=0)
        self.log.debug("self.kpoints.shape (after concat.): " +
                       str(self.kpoints.shape))
        self.sx = np.concatenate(rotations[1], axis=0)
        self.sy = np.concatenate(rotations[2], axis=0)
        self.sz = np.concatenate(rotations[3], axis=0)
        # the bands and proj. character also need to be enlarged
        bandsChar = [(self.bands, self.character) for i in range(order)]
        bandsChar = zip(*bandsChar)
        self.bands = np.concatenate(bandsChar[0], axis=0)
        self.character = np.concatenate(bandsChar[1], axis=0)
        self.log.debug("RotSymmZ:...Done")

        return

    def Translate(self, newOrigin):
        """Centers the Kpoints at newOrigin, newOrigin is either and index (of
       some Kpoint) or the cartesian coordinates of one point in the
       reciprocal space.

        """
        self.log.debug("Translate():  ...")
        if len(newOrigin) == 1:
            newOrigin = int(newOrigin[0])
            newOrigin = self.kpoints[newOrigin]
        # Make sure newOrigin is a numpy array
        newOrigin = np.array(newOrigin)
        self.log.debug("newOrigin: " + str(newOrigin))
        self.kpoints -= newOrigin
        self.log.debug("new Kpoints:\n" + str(self.kpoints))
        self.log.debug("Translate(): ...Done")
        return


def scriptCat(args):
    if args.quiet is False:
        print("Concatenating:")
        print("Input         : ", ', '.join(args.inFiles))
        print("Output        : ", args.outFile)
        if args.gz:
            print("out compressed: True")
        if args.verbose > 2:
            args.verbose = 2
        if args.verbose:
            print("verbosity     : ", args.verbose)

    if args.gz and args.outFile[-3:] is not '.gz':
        args.outFile += '.gz'
        if args.quiet is False:
            print(".gz extension appended to the outFile")

    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]
    handler = UtilsProcar(loglevel)
    handler.MergeFiles(args.inFiles, args.outFile, gzipOut=args.gz)
    return


def scriptFilter(args):
    if args.quiet is False:
        print("Input file  :", args.inFile)
        print("Output file :", args.outFile)
    if args.verbose:
        if args.verbose > 2:
            args.verbose = 2
        print("atoms       :", args.atoms)
        if args.atoms:
            print("human     :", args.human)
        print("orbitals  :", args.orbitals)
        if args.orbitals:
            print("orb. names  :", args.orbital_names)
        print("bands       :", args.bands)
        print("spins       :", args.spin)

    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]
    FileFilter = ProcarFileFilter(args.inFile, args.outFile, loglevel=loglevel)

    if args.atoms:
        if args.quiet is False:
            print("Manipulating the atoms")

        if args.human:
            args.atoms = [[y - 1 for y in x] for x in args.atoms]
            if args.verbose:
                print("new atoms list :", args.atoms)

        # Now just left to call the driver member
        FileFilter.FilterAtoms(args.atoms)

    elif args.orbitals:
        if args.quiet is False:
            print("Manipulating the orbitals")
        # If orbitals orbital_names is None, it needs to be filled
        if args.orbital_names is None:
            args.orbital_names = ["o" + str(x) for x in range(len(args.orbitals))]
            if args.quiet is False:
                print("New orbitals names (default): ", args.orbital_names)
        # testing if makes sense
        if len(args.orbitals) != len(args.orbital_names):
            raise RuntimeError("length of orbitals and orbitals names do not match")

        FileFilter.FilterOrbitals(args.orbitals, args.orbital_names)

    elif args.bands:
        if args.quiet is False:
            print("Manipulating the bands")

        bmin = args.bands[0]
        bmax = args.bands[1]
        if bmax < bmin:
            bmax, bmin = bmin, bmax
            if args.quiet is False:
                print("New bands limits: ", bmin, " to ", bmax)

        FileFilter.FilterBands(bmin, bmax)

    elif args.spin:
        if args.quiet is False:
            print("Manipulating the spin")

        FileFilter.FilterSpin(args.spin)

    return


def scriptBandsplot(args):
    # First handling the options, to get feedback to the user and check
    # that the input makes sense.
    # It is quite long
    if args.atoms is None:
        args.atoms = [-1]
        if args.human is True:
            print("WARNING: `--human` option given without atoms list!")
            print("--human will be set to False (ignored)\n ")
            args.human = False
    if args.orbitals is None:
        args.orbitals = [-1]

    if args.verbose:
        print("Script initiated")
        print("input file    : ", args.file)
        print("Mode          : ", args.mode)

        print("spin comp.    : ", args.spin)
        print("atoms list.   : ", args.atoms)
        print("orbs. list.   : ", args.orbitals)

    if args.fermi is None and args.outcar is None:
        print("WARNING: Fermi Energy not set! ")
        print("You should use '-f' or '--outcar'")
        print("The zero of energy is arbitrary\n")
        args.fermi = 0

    if args.verbose:
        print("Fermi Ener.   : ", args.fermi)
        print("Energy range  : ", args.elimit)

        if args.mask is not None:
            print("masking thres.: ", args.mask)

        print("Colormap      : ", args.cmap)
        print("MarkerSize    : ", args.markersize)

        print("Permissive    : ", args.permissive)
        if args.permissive and args.quiet is False:
            print("INFO: Permissive flag is on! Be careful")
        print("vmax          : ", args.vmax)
        print("vmin          : ", args.vmin)
        print("grid enabled  : ", args.grid)
        if args.human is not None:
            print("human         : ", args.human)
        print("Savefig       : ", args.savefig)
        print("kticks        : ", args.kticks)
        print("knames        : ", args.knames)
        print("title         : ", args.title)

        print("outcar        : ", args.outcar)

    # If ticks and names are given we should use them#
    if args.kticks is not None and args.knames is not None:
        ticks = zip(args.kticks, args.knames)
    elif args.kticks is not None:
        ticks = zip(args.kticks, args.kticks)
    else:
        ticks = None

    # The spin argument should be a number (index of an array), or
    # 'st'. In the last case it will be handled separately (later)
    args.spin = {'0': 0, '1': 1, '2': 2, '3': 3, 'st': 'st'}[args.spin]

    if args.verbose > 2:
        args.verbose = 2
    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]

    # The second part of this function is parse/select/use the data in
    # OUTCAR (if given) and PROCAR

    # first parse the outcar if given, to get Efermi and Reciprocal lattice
    recLat = None
    if args.outcar:
        outcarparser = UtilsProcar(loglevel=loglevel)
        if args.fermi is None:
            args.fermi = outcarparser.FermiOutcar(args.outcar)
            if args.quiet is False:
                print("INFO: Fermi energy found in outcar file = " + str(args.fermi))
        recLat = outcarparser.RecLatOutcar(args.outcar)

    # parsing the PROCAR file
    procarFile = ProcarParser(loglevel=loglevel)
    procarFile.readFile(args.file, args.permissive, recLat)

    # processing the data, getting an instance of the class that reduces the data
    data = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)

    # handling the spin, `args.spin='st'` is not straightforward, needs
    # to calculate the k vector and its normal. Other `args.spin` values
    # are trivial.
    if args.spin is 'st':
        # two `ProcarSelect` instances, to store temporal values: spin_x, spin_y
        dataX = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)
        dataX.selectIspin([1])
        dataX.selectAtoms(args.atoms, fortran=args.human)
        dataX.selectOrbital(args.orbitals)
        dataY = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)
        dataY.selectIspin([2])
        dataY.selectAtoms(args.atoms, fortran=args.human)
        dataY.selectOrbital(args.orbitals)
        # getting the signed angle of each K-vector
        angle = np.arctan2(dataX.kpoints[:, 1], (dataX.kpoints[:, 0] + 0.000000001))
        sin = np.sin(angle)
        cos = np.cos(angle)
        sin.shape = (sin.shape[0], 1)
        cos.shape = (cos.shape[0], 1)
        # #print sin, cos
        # storing the spin projection into the original array
        data.spd = -sin * dataX.spd + cos * dataY.spd
    else:
        data.selectIspin([args.spin])
        data.selectAtoms(args.atoms, fortran=args.human)
        data.selectOrbital(args.orbitals)

    # Plotting the data
    assert (data.bands is not None)
    data.bands = (data.bands.transpose() - np.array(args.fermi)).transpose()
    plot = ProcarPlot(data.bands, data.spd, data.kpoints)

    # start of mode dependent options
    if args.mode == "scatter":
        plot.scatterPlot(mask=args.mask, size=args.markersize,
                         cmap=args.cmap, vmin=args.vmin,
                         vmax=args.vmax, marker=args.marker, ticks=ticks)

        plt.ylabel(r"Energy [eV]")
        if args.elimit is not None:
            plt.ylim(args.elimit)

    elif args.mode == "plain":
        plot.plotBands(args.markersize, marker=args.marker, ticks=ticks)
        plt.ylabel(r"Energy [eV]")
        if args.elimit:
            plt.ylim(args.elimit)

    elif args.mode == "parametric":
        plot.parametricPlot(cmap=args.cmap, vmin=args.vmin, vmax=args.vmax,
                            ticks=ticks)
        plt.ylabel(r"Energy [eV]")
        if args.elimit is not None:
            plt.ylim(args.elimit)

    elif args.mode == "atomic":
        plot.atomicPlot(cmap=args.cmap, vmin=args.vmin, vmax=args.vmax)
        plt.ylabel(r"Energy [eV]")
        if args.elimit is not None:
            plt.ylim(args.elimit)
    # end of mode dependent options

    if args.grid:
        plt.grid()

    if args.title:
        plt.title(args.title)

    if args.savefig:
        plt.savefig(args.savefig, bbox_inches=0)
    else:
        plt.show()

    return


def scriptFermi2D(args):
    if args.atoms is None:
        args.atoms = [-1]
        if args.human is True:
            print("WARNING: `--human` option given without atoms list!!!!!")

    if args.orbitals is None:
        args.orbitals = [-1]

    if args.rec_basis is not None:
        args.rec_basis = np.array(args.rec_basis)
        args.rec_basis.shape = (3, 3)

    if len(args.translate) != 3 and len(args.translate) != 1:
        print("Error: --translate option is invalid! (", args.translate, ")")
        raise RuntimeError("invalid option --translate")

    if args.quiet is False:
        print("file            : ", args.file)
        print("atoms           : ", args.atoms)
        print("orbitals        : ", args.orbitals)
        print("spin comp.      : ", args.spin)
        print("energy          : ", args.energy)
        print("fermi energy    : ", args.fermi)
        print("Rec. basis      : ", args.rec_basis)
        print("rot. symmetry   : ", args.rot_symm)
        print("origin (trasl.) : ", args.translate)
        print("rotation        : ", args.rotation)
        print("masking thres.  : ", args.mask)
        print("save figure     : ", args.savefig)
        print("outcar          : ", args.outcar)
        print("st              : ", args.st)
        print("no_arrows       : ", args.noarrow)
    if args.verbose > 2:
        args.verbose = 2
    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]

    # first parse the outcar, if given
    if args.rec_basis is None and args.outcar:
        outcarparser = UtilsProcar(loglevel=loglevel)
        if args.fermi is None:
            args.fermi = outcarparser.FermiOutcar(args.outcar)
            if args.quiet is False:
                print("Fermi energy found in outcar file = " + str(args.fermi))
        args.rec_basis = outcarparser.RecLatOutcar(args.outcar)
    # Reciprocal lattices are needed!
    elif args.rec_basis is None and args.outcar is None:
        print("ERORR: Reciprocal Lattice is needed, use --rec_basis or --outcar")
        raise RuntimeError("Reciprocal Lattice not found")

    # parsing the file
    procarFile = ProcarParser(loglevel)
    # permissive incompatible with Fermi surfaces
    procarFile.readFile(args.file, permissive=False, recLattice=args.rec_basis)

    if args.st is not True:
        # processing the data
        data = ProcarSelect(procarFile, loglevel=loglevel)
        data.selectIspin([args.spin])
        # fortran flag is equivalent to human,
        # but the later seems more human-friendly
        data.selectAtoms(args.atoms, fortran=args.human)
        data.selectOrbital(args.orbitals)
    else:
        # first get the sdp reduced array for all spin components.
        stData = []
        for i in [1, 2, 3]:
            data = ProcarSelect(procarFile, loglevel=loglevel)
            data.selectIspin([i])
            data.selectAtoms(args.atoms, fortran=args.human)
            data.selectOrbital(args.orbitals)
            stData.append(data.spd)

    # Once the PROCAR is parsed and reduced to 2x2 arrays, we can apply
    # symmetry operations to unfold the Brillouin Zone
    kpoints = data.kpoints
    bands = data.bands
    character = data.spd
    if args.st is True:
        sx, sy, sz = stData[0], stData[1], stData[2]
        symm = ProcarSymmetry(kpoints, bands, sx=sx, sy=sy, sz=sz,
                              loglevel=logging.DEBUG)
    else:
        symm = ProcarSymmetry(kpoints, bands, character=character,
                              loglevel=logging.DEBUG)

    symm.Translate(args.translate)
    symm.GeneralRotation(args.rotation[0], args.rotation[1:])
    symm.RotSymmetryZ(args.rot_symm)

    # visual the data
    print("Bands will be shifted by the Fermi energy = ", args.fermi)
    fs = FermiSurface(symm.kpoints, symm.bands - args.fermi, character,
                      loglevel=loglevel)
    fs.FindEnergy(args.energy)

    if not args.st:
        fs.Plot(mask=args.mask, interpolation=300)
    else:
        fs.st(sx=symm.sx, sy=symm.sy, sz=symm.sz, noarrow=args.noarrow, spin=args.spin)

    if args.savefig:
        plt.savefig(args.savefig)
    else:
        plt.show()

    return


def scriptVector(args):
    if args.quiet is False:
        print("Input File    : ", args.infile)
        print("Bands         : ", args.bands)
        print("Energy        : ", args.energy)
        print("Fermi         : ", args.fermi)
        print("outcar        : ", args.outcar)
        print("atoms         : ", args.atoms)
        print("orbitals      : ", args.orbitals)
        print("scale factor  : ", args.scale)

    if args.bands is [] and args.energy is None:
        raise RuntimeError("You must provide the bands or energy.")
    if args.fermi is None and args.outcar is None:
        print("WARNING: Fermi's Energy not set")
    if args.verbose > 2:
        args.verbose = 2
    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]

    # first parse the outcar if given
    recLat = None  # Will contain reciprocal vectors, if necessary
    if args.outcar:
        outcarparser = UtilsProcar(loglevel=loglevel)
        if args.fermi is None:
            args.fermi = outcarparser.FermiOutcar(args.outcar)
            if args.quiet is False:
                print("Fermi energy found in outcar file = " + str(args.fermi))
        recLat = outcarparser.RecLatOutcar(args.outcar)

    if args.atoms is None:
        args.atoms = [-1]
    if args.orbitals is None:
        args.orbitals = [-1]

    # parsing the file
    procarFile = ProcarParser(loglevel=loglevel)
    procarFile.readFile(args.infile, recLattice=recLat)

    # processing the data
    sx = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)
    sx.selectIspin([1])
    sx.selectAtoms(args.atoms)
    sx.selectOrbital(args.orbitals)

    sy = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)
    sy.selectIspin([2])
    sy.selectAtoms(args.atoms)
    sy.selectOrbital(args.orbitals)

    sz = ProcarSelect(procarFile, deepCopy=True, loglevel=loglevel)
    sz.selectIspin([3])
    sz.selectAtoms(args.atoms)
    sz.selectOrbital(args.orbitals)

    x = sx.kpoints[:, 0]
    y = sx.kpoints[:, 1]
    z = sx.kpoints[:, 2]

    # if energy was given I need to find the bands indexes crossing it
    if args.energy is not None:
        FerSurf = FermiSurface(sx.kpoints, sx.bands - args.fermi, sx.spd, loglevel)
        FerSurf.FindEnergy(args.energy)
        args.bands = list(FerSurf.useful[0])
        print("Bands indexes crossing Energy  ", args.energy, ", are: ", args.bands)

    from mayavi import mlab

    mlab.figure(bgcolor=(1, 1, 1))

    for band in args.bands:
        # z = sx.bands[:,band]-args.fermi
        u = sx.spd[:, band]
        v = sy.spd[:, band]
        w = sz.spd[:, band]
        scalar = w

        vect = mlab.quiver3d(x, y, z, u, v, w, scale_factor=args.scale,
                             scale_mode='vector', scalars=scalar, mode='arrow',
                             colormap='jet')
        vect.glyph.color_mode = 'color_by_scalar'
        vect.scene.parallel_projection = True
        vect.scene.z_plus_view()

        # tube= mlab.plot3d(x,y,z, tube_radius=0.0050, color=(0.5,0.5,0.5))
    mlab.show()


def script_repair(args):
    if args.quiet is False:
        print("Input File    : ", args.infile)
        print("Output File   : ", args.outfile)

    if args.verbose > 2:
        args.verbose = 2
    loglevel = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}[args.verbose]

    # parsing the file
    handler = UtilsProcar(loglevel=loglevel)
    handler.ProcarRepair(args.infile, args.outfile)


if __name__ == "__main__":
    import argparse
    from argparse import RawTextHelpFormatter

    description = ("procar.py: a python library/script to manipulate and "
                   "plot VASP's PROCAR files.")
    # A top-level parser

    parser = argparse.ArgumentParser(description=description)
    subparsers = parser.add_subparsers(help='sub-command')

    # concatenation
    phelp = ("concatenation of PROCARs files, they should be compatible (ie: "
             "joining parts of a large bandstructure calculation).")
    parserCat = subparsers.add_parser('cat', help=phelp)

    phelp = "Input files. They can be compressed"
    parserCat.add_argument('inFiles', nargs='+', help=phelp)

    phelp = "Output file."
    parserCat.add_argument('outFile', help=phelp)

    phelp = ("Writes a gzipped outfile (if needed a .gz extension automatically "
             "will be added)")
    parserCat.add_argument('--gz', help=phelp, action='store_true')

    VerbCat = parserCat.add_mutually_exclusive_group()
    VerbCat.add_argument("-v", "--verbose", action="count", default=0)
    VerbCat.add_argument("-q", "--quiet", action="store_true")

    parserCat.set_defaults(func=scriptCat)

    # filter
    phelp = ("Filters (manipulates) the data of the input file (PROCAR-like) and"
             " it yields a new file (PROCAR-like too) with the changes. This "
             "method can do only one manipulation at time (ie: spin, atoms, "
             "bands or orbitals).")
    parserFilter = subparsers.add_parser('filter', help=phelp)

    phelp = "Input file. Can be compressed"
    parserFilter.add_argument('inFile', help=phelp)

    phelp = 'Output file.'
    parserFilter.add_argument('outFile', help=phelp)

    VerbFilter = parserFilter.add_mutually_exclusive_group()
    VerbFilter.add_argument("-v", "--verbose", action="count", default=0)
    VerbFilter.add_argument("-q", "--quiet", action="store_true")

    OptFilter = parserFilter.add_mutually_exclusive_group()
    phelp = ("List of atoms to group (add) as a new single entry. Each group of"
             " atoms should be specified in a different `--atoms` option. "
             "Example: `procar.py filter in out -a 0 1 -a 2` will group the 1st"
             " and 2nd atoms, while keeping the 3rd atom in `out` (any atom "
             "beyond the 3rd will be discarded). Mind the last atomic field "
             "present on a PROCAR file, is not an atom, is the 'tot' value (sum"
             " of all atoms), this field always is included in the outfile and "
             "it always is the 'tot' value from infile, regardless the selection"
             " of atoms.")
    OptFilter.add_argument("-a", "--atoms", type=int, nargs='+', action='append',
                           help=phelp)

    phelp = ("List of orbitals to group as a single entry. Each group of "
             "orbitals needs a different `--orbitals` list. By instance, to "
             "group orbitals in 's','p', 'd' it is needed `-o 0 -o 1 2 3 -o 4 5 "
             "6 7 8`. Where 0=s, 1,2,3=px,py,pz, 4...9=dxx...dyz. Mind the last "
             "value (aka) 'tot' always is written.")
    OptFilter.add_argument("-o", "--orbitals", help=phelp, type=int, nargs='+',
                           action='append')

    phelp = ("Keeps only the bands between `min` and `max` indexes. To keep the "
             "bands from 120 to 150 you should give `-b 120 150 `. It is not "
             "obvious which indexes are in the interest region, therefore I "
             "recommend you trial and error ")
    OptFilter.add_argument("-b", "--bands", help=phelp, type=int, nargs=2)

    phelp = ("Which spin components should be written: 0=density, 1,2,3=Sx,Sy,Sz."
             " They are not averaged.")
    OptFilter.add_argument("-s", "--spin", help=phelp, type=int, nargs='+')

    phelp = ("enable to give atoms list in a more human, 1-based order (say the"
             " 1st is 1, 2nd is 2 and so on ). Mind: this only holds for atoms.")
    parserFilter.add_argument("--human", help=phelp, action="store_true")

    phelp = ("List of names of new 'orbitals' to appear in the new file, eg. "
             "(`--orbital_names s p d` for a 's', 'p', 'd'). Only meaningful "
             "when manipulating the orbitals, ie: using `-o` ")
    parserFilter.add_argument("--orbital_names", help=phelp, nargs='+')

    parserFilter.set_defaults(func=scriptFilter)

    # bandsplot
    phelp = ("Bandstructure plot. This kind of plot can be fairly complex, "
             "therefore its worth to explore all options. If the file is large "
             "(>100MB) you should consider use the 'procar.py filter' command "
             "before. Mind not all option are meaningful (even considered) for "
             "all '--modes'")
    parserBandsplot = subparsers.add_parser('bandsplot', help=phelp,
                                            formatter_class=RawTextHelpFormatter)

    phelp = "Input file. It can be compressed"
    parserBandsplot.add_argument('file', help=phelp)

    VerbBandsPlot = parserBandsplot.add_mutually_exclusive_group()
    VerbBandsPlot.add_argument("-v", "--verbose", action="count", default=0)
    VerbBandsPlot.add_argument("-q", "--quiet", action="store_true")

    phelp = ("Mode of plot, how do you like to plot the bands? Available modes:\n"
             "-m  scatter : is a points plot with the color given by the chosen\n"
             "  projection. It produces a rather heavy pdf file.\n\n"
             "-m  parametric : like scatter, but with lines instead of points \n"
             "  (bands crossings are not handled, and some  unphysical 'jumps' \n"
             " can be present). Sligthy smaller PDF size.\n\n"
             "-m plain : is a featureless bandstructure ignoring all data about\n"
             "  projections. Rather light-weight\n\n"
             "-m atomic : For non-periodic system, like molecules, rather ugly \n"
             "  but useful to visualize energy level. Only 1 K-point!\n\n")
    choices = ["scatter", "plain", "parametric", "atomic"]
    parserBandsplot.add_argument("-m", "--mode", help=phelp, default="scatter",
                                 choices=choices)

    phelp = ("Spin component to be used (default -s 0): \n\n"
             "Non-polarized calculations density is '-s 0', just ignore it.\n\n"
             "Spin-Polarized (collinear) calculation: \n"
             "-s 0 are the unpolarized bands, the spin-polarization is ignored.\n"
             "-s 1 'spin-polarized' bands, the character of 'up' bands positive,\n"
             "  but negative for 'down' bands, this means that the color of \n"
             "  'down' is negative. Useful together with '--cmap seismic'.\n\n"
             "Non-collinear calculation: \n"
             "-s 0 : density, ie: Spin-orbit-coupling but don't care of spin.\n"
             "-s 1 : Sx, projection along 'x' quantization axis, see SAXIS flag\n"
             "  in the VASP manual\n"
             "-s 2 : Sy, projection along 'y' quantization axis.\n"
             "-s 3 : Sy, projection along 'z' quantization axis.\n"
             "-s st : Spin-texture perpendicular in the plane (kx,ky) to each\n "
             "(kx,ky) vector. Useful for Rashba-like states in surfaces. Use\n "
             "'--cmap seismic'\n\n ")
    parserBandsplot.add_argument("-s", "--spin", choices=['0', '1', '2', '3', 'st'],
                                 default='0', help=phelp)

    phelp = ("List of rows (atoms) to be used. This list refers to the rows of\n"
             "(each block of) your PROCAR file. If you haven't manipulated your\n"
             "PROCAR (eg: with the '-a' option of 'filter' mode) each row\n"
             "correspond to the respective atom in the POSCAR.\n\n"
             "Mind: This list is 0-based, ie: the 1st atom is 0, the 2nd is 1,\n"
             "  and so on. If you need to be treated like a human, specify '-u'\n"
             "or '--human' and 1st->1, 2nd->2, etc.\n\n"
             "Example:\n"
             "-a 0 2 :  select the 1st  and 3rd. rows (likely 1st and 3rd atoms)"
             "\n\n")
    parserBandsplot.add_argument("-a", "--atoms", type=int, nargs='+',
                                 help=phelp)

    phelp = ("Orbitals index(es) to be used, take a look to the PROCAR file, \n"
             "they are 's py pz px ...', then s->0, py->1, pz->2 and so on. \n"
             "Note that indexes begin at 0!. Its default is the last field (ie:\n"
             "'tot', did you saw the PROCAR?). Some examples:\n\n"
             "-o 0 : s-orbital (unless you modified the orbitals, eg. 'filter')\n"
             "-o 1 2 3 : py+pz+px (unless you modified the orbitals)\n"
             "-o 4 5 6 7 8 : all the d-orbitasl (unless...)\n"
             "-o 2 6 : pz+dzz (did you look at the PROCAR?)\n\n ")
    parserBandsplot.add_argument("-o", "--orbitals", type=int, nargs='+',
                                 help=phelp)

    phelp = ("Set the Fermi energy (or any reference energy) as the zero energy.\n"
             "See '--outcar', avoids to give it explicitly. A list of \n"
             "k-dependant 'fermi-like energies' are also accepted (useful to\n"
             "compare different systems in one PROCAR made by hand). \n\n"
             "Mind: The Fermi energy MUST be the one from the self-consistent\n"
             "calculation, not from a Bandstructure calculation!\n\n")
    parserBandsplot.add_argument("-f", "--fermi", type=float, help=phelp,
                                 nargs='+')

    phelp = ("Min/Max energy to be ploted. Example:\n "
             "--elimit -1 1 : From -1 to 1 around Fermi energy (if given)\n\n")
    parserBandsplot.add_argument("--elimit", type=float, nargs=2, help=phelp)

    phelp = ("If given, it masks(hides) bands with values lowers than 'mask'.\n"
             "It is useful to remove 'unwanted' bands. For instance, if you\n"
             "project the bandstructure on a 'surface' atom -with the default\n"
             "colormap- some white points can appear, they are bands with \n"
             "almost no contribution to the 'surface': no physics but they \n"
             "still look ugly, to hide those bands use '--mask 0.08' (or some \n"
             "other small value). Mind: it works with the absolute value of\n"
             "projection (no problem with spin polarization)\n\n")
    parserBandsplot.add_argument("--mask", type=float, help=phelp)

    phelp = ("Size of markers, if used. Each mode has it own scale,\n"
             "just test them\n\n")
    parserBandsplot.add_argument("--markersize", type=float, help=phelp,
                                 default=10)

    phelp = ("Do you find ugly the colors I choose? Depending on the context\n"
             "another color scheme can be more meaningful. It is super-easy to\n"
             "change just specify the name of the color map after '--cmap'\n"
             "Don't you know the names? I don't know either!, but GOOGLE knows\n"
             "them, just search 'python colormap'). For instance:\n\n"
             "--cmap  seismic : blue->white->red, useful to see the \n"
             "  spin-polarization of a band (it will blueish or reddish)\n"
             "  depending of spin channel\n"
             "--cmap  seismic_r : the 'seismic' colormap, but reversed.\n\n")
    parserBandsplot.add_argument("--cmap", help=phelp, default="hot_r")

    phelp = ("Do you want to Normalize the plots to the same scale of colors\n"
             "(ie: the numbers on the bar at the right), just try '--vmax'\n\n"
             "--vmax 1 : If you are looking for the s, p or d character.\n"
             "--vmax 0.2: If you want to capture some tiny effect, eg: s-band of\n"
             "  a impurity on a metal\n\n")
    parserBandsplot.add_argument("--vmax", type=float, help=phelp)

    phelp = ("Like '--vmax' (see '--vmax'), However, for spin-polarized \n"
             "(collinear or not) you can set it to a negative value. Actually\n"
             "you can do it for a non-spin-polarized calculation and the \n"
             "effect will be a 'stretching' of the color scheme, try it.\n\n")
    parserBandsplot.add_argument("--vmin", type=float, help=phelp)

    phelp = "switch on/off the grid. Default is 'on'\n\n"
    parserBandsplot.add_argument("--grid", type=bool, help=phelp, default=True)

    phelp = ("set the marker shape, ie: 'o'=circle, 's'=square,\ '-'=line\n"
             "(only mode `plain`, other symbols: google pyplot markers)\n\n")
    parserBandsplot.add_argument("--marker", type=str, help=phelp, default='o')

    phelp = ("Some fault tolerance for ill-formatted files (stupid fortran)\n"
             "But be careful, something could be messed up and don't work (at\n"
             "least as expected). Length of K-points paths will be ignored\n\n")
    parserBandsplot.add_argument("--permissive", help=phelp, action='store_true')

    phelp = ("Enable human-like 1-based order (ie 1st is 1, 2nd is 2, and so\n"
             "on). Mind: this only works for atoms, not for orbitals or spin\n\n")
    parserBandsplot.add_argument("-u", "--human", help=phelp, action="store_true")

    phelp = ("Saves the figure, instead of display it on screen. Anyway, you can\n"
             "save from the screen too. Any file extension supported by\n "
             "`matplotlib.savefig` is valid (if you are too lazy to google it,\n"
             "trial and error also works fine)\n\n")
    parserBandsplot.add_argument('--savefig', help=phelp)

    phelp = ("list of ticks along the kpoints axis (x axis). For instance a\n"
             "bandstructure G-X-M with 10 point by segment should be:\n "
             "--kticks 0 9 19\n\n")
    parserBandsplot.add_argument('--kticks', help=phelp, nargs='+', type=int)

    phelp = ("Names of the points given in `--kticks`. In the `kticks` example\n"
             "they should be `--knames \"\$Gamma\$\" X M`. As you can see \n"
             "LaTeX stuff works with a minimal mess (extra \\s)\n\n")
    parserBandsplot.add_argument('--knames', help=phelp, nargs='+',
                                 type=str)

    phelp = ("Title, to use several words, use quotation marks\"\" or ''. Latex\n"
             " works if you scape the special characteres, ie: $\\alpha$ -> \n"
             "\$\\\\alpha\$")
    parserBandsplot.add_argument('-t', '--title', help=phelp, type=str)

    phelp = ("OUTCAR file where to find the reciprocal lattice vectors and\n "
             "perhaps E_fermi.\n"
             "Mind: '--fermi' has precedence, remember that the E-fermi should\n"
             "correspond to a self-consistent run, not to a bandstructure!\n"
             "(however, the basis and reciprocal vectors will be safe from the\n"
             "non-self-consistentness)\n\n")
    parserBandsplot.add_argument('--outcar', help=phelp)

    parserBandsplot.set_defaults(func=scriptBandsplot)

    # Fermi
    phelp = ("Plot the Fermi surface for a 2D Brillouin zone (layer-wise)"
             " along z and just perpendicular to z! (actually k_z)")
    parserFermi2D = subparsers.add_parser('fermi2D', help=phelp)

    phelp = "Input file. It Can be compressed"
    parserFermi2D.add_argument('file', help=phelp)

    VerbFermi2D = parserFermi2D.add_mutually_exclusive_group()
    VerbFermi2D.add_argument("-v", "--verbose", action="count", default=0)
    VerbFermi2D.add_argument("-q", "--quiet", action="store_true")

    phelp = ("Spin component to be used: for non-polarized calculations density "
             "is '-s 0'. For spin polarized case test '-s 0' (ignore the spin) or"
             " '-s 1' (assing a sign to the spin channel). For "
             "non-collinear stuff you can use '-s 0', '-s 1', '-s 2', -s "
             "3 for the magnitude, x, y, z components of your spin "
             "vector, in this case you really want to see spin textures "
             "'--st'. Default: s=0")
    parserFermi2D.add_argument("-s", "--spin", type=int, choices=[0, 1, 2, 3],
                               default=0, help=phelp)

    phelp = ("List of atoms to be used (0-based): ie. '-a 0 2' to select the 1st"
             " and 3rd. It defaults to the last one (-a -1 the 'tot' entry)")
    parserFermi2D.add_argument("-a", "--atoms", type=int, nargs='+', help=phelp)

    phelp = ("Orbital index(es) to be used 0-based. Take a look to the PROCAR "
             "file. Its default is the last field (ie: 'tot'). From a standard "
             "PROCAR: `-o 0`='s', `-o 1 2 3`='p', `-o 4 5 6 7 8`='d'.")
    parserFermi2D.add_argument("-o", "--orbitals", type=int, nargs='+',
                               help=phelp)

    phelp = ("Energy for the surface. To plot the Fermi surface at Fermi Energy "
             "`-e 0`")
    parserFermi2D.add_argument("-e", "--energy", help=phelp, type=float,
                               required=True)

    phelp = ("Set the Fermi energy (or any reference energy) as zero. To get it "
             "you should `grep E-fermi` the self-consistent OUTCAR. See "
             "`--outcar`")
    parserFermi2D.add_argument("-f", "--fermi", help=phelp, type=float)

    phelp = ("reciprocal space basis vectors. 9 number are required b1x b1y ... "
             " b3z. This option is quite involved, so I recommend you to use "
             "`--outcar`")
    parserFermi2D.add_argument("--rec_basis", help=phelp, type=float, nargs=9)

    phelp = ("Apply a rotational symmetry to unfold the Kpoints found. If your"
             "PROCAR only has a portion of the Brillouin Zone, you may want to "
             "plot the FULL BZ (ie: a Dirac cone at Gamma will look like a cone "
             "and not like a segment of circle). Supported rotations are "
             "1,2,3,4,6. All of them along Z and centered at Gamma. Consider to "
             "'--translate' your cell to rotate with other origin. This is the "
             "last symmetry operation to be performed.")
    parserFermi2D.add_argument("--rot_symm", help=phelp, type=int, default=1)

    phelp = ("Translate your mesh to the specified point. The point can be 3 "
             "coordinates (numbers) or the index of one K-point (zero-based, as "
             "usual). This is the first symmetry operation to be performed "
             "(i.e. rotations will take this point as the origin).")
    parserFermi2D.add_argument("--translate", help=phelp, nargs='+',
                               default=[0, 0, 0])

    phelp = ("A general rotation is applied to the data in the PROCAR. While this "
             " script has a large bias to work on the 'xy' plane, with this option"
             " you can rotate your whole PROCAR to fit the 'xy' plane. A rotation "
             "is composed by one angle plus one fixed axis, eg: '--rotation 90 1 0"
             " 0' is 90 degrees along the x axis, this changes the 'xy'->'xz'. The"
             " rotation is performed after the translation and before applying "
             "rot_symm. ")
    parserFermi2D.add_argument("--rotation", help=phelp, type=float,
                               nargs=4, default=[0, 0, 0, 1])

    phelp = ("enable to give atoms list"
             " in a more human, 1-based way (say the 1st is 1,"
             " 2nd is 2 and so on )")
    parserFermi2D.add_argument("-u", "--human", help=phelp, action="store_true")

    phelp = ("If set, masks(hides)"
             " values lowers than it. Useful to remove "
             "unwanted bands.")
    parserFermi2D.add_argument("--mask", type=float, help=phelp)

    phelp = ("Saves the figure, instead of "
             "display it on screen. Anyway, you can save from"
             " the screen too. Any file extension supported by"
             "matplotlib.savefig is valid (if you are too lazy"
             " to google it, trial and error also works)")
    parserFermi2D.add_argument('--savefig', help=phelp)

    phelp = ("OUTCAR file where to find the reciprocal lattice vectors and "
             "perhaps E_fermi. Mind: '--fermi' has precedence, remember that the"
             " E-fermi should correspond to a self-consistent run, not a "
             "bandstructure! (however, this is irrelevant for basis vectors)")
    parserFermi2D.add_argument('--outcar', help=phelp)

    phelp = ("Plot of the spin texture (ie: spin arrows) on the Fermi's surface."
             " This option works quite indepentent of another options.")
    parserFermi2D.add_argument('--st', help=phelp, action='store_true')

    phelp = ("Plot of the spin texture without arrows (just intensity) for a "
             "given spin direction on the Fermi's surface.  This option works"
             "quite indepentent of another options but needs to set '--st' and "
             "'--spin'.")
    parserFermi2D.add_argument('--noarrow', help=phelp, action='store_true')

    parserFermi2D.set_defaults(func=scriptFermi2D)

    # vector
    # only for non-collinear bandstructures

    phelp = "Plots the bands (some of them) as a vector field. It uses Mayavi."
    parserVector = subparsers.add_parser('Vector', help=phelp)

    phelp = "Input file. Can be compressed"
    parserVector.add_argument('infile', help=phelp)

    VerbVector = parserVector.add_mutually_exclusive_group()
    VerbVector.add_argument("-v", "--verbose", action="count", default=0)
    VerbVector.add_argument("-q", "--quiet", action="store_true")

    EnVector = parserVector.add_mutually_exclusive_group()
    phelp = "Band(s) to plot, by bands index"
    EnVector.add_argument("-b", "--bands", nargs='+', type=int, help=phelp)
    phelp = "Band(s) to plot, by Energy"
    EnVector.add_argument("-e", "--energy", type=float, help=phelp)

    phelp = "Fermi level. All energies will be referred to it as zero"
    parserVector.add_argument('-f ', '--fermi', type=float, help=phelp)

    phelp = "List of atoms to consider (add their contribution)"
    parserVector.add_argument('-a ', '--atoms', type=int, nargs='+', help=phelp)

    phelp = "List of orbital to consider (add their proj.)"
    parserVector.add_argument('-o ', '--orbitals', type=int, nargs='+',
                              help=phelp)

    phelp = ("OUTCAR file where to grep the reciprocal lattice vectors and "
             "perhaps E_fermi. Mind: '--fermi' has precedence, remember that "
             "the E-fermi should correspond to a self-consistent run, not a "
             "bandstructure!")
    parserVector.add_argument('--outcar', help=phelp)

    phelp = "Scale factor, to avoid too large arrows"
    parserVector.add_argument('--scale', type=float, default=0.1, help=phelp)

    parserVector.set_defaults(func=scriptVector)

    # repair
    phelp = "Attemp to repair a Procar file form some fixed format problems."
    parserRepair = subparsers.add_parser('repair', help=phelp)

    phelp = "Input file. Can be compressed"
    parserRepair.add_argument('infile', help=phelp)

    phelp = "Output file"
    parserRepair.add_argument('outfile', help=phelp)

    VerbRepair = parserRepair.add_mutually_exclusive_group()
    VerbRepair.add_argument("-v", "--verbose", action="count", default=0)
    VerbRepair.add_argument("-q", "--quiet", action="store_true")

    parserRepair.set_defaults(func=script_repair)

    args = parser.parse_args()
    args.func(args)
