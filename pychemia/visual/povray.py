import numpy as np
from pychemia.utils.periodic import atomic_number, covalent_radius, cpk_colors


class StructurePovray:
    def __init__(self, structure):

        self.structure = structure
        self.distance = 10

    def create_pov(self):

        ret = """
#version 3.7;
#include "colors.inc"    // The include files contain
#include "stones.inc"    // pre-defined scene elements
#include "glass.inc"
background{rgb 0}

"""
        if self.structure.is_crystal:
            self.distance = max(self.structure.lattice.lengths)
        else:
            self.distance = 10

        ret += "#declare r=%7.3f;\n #declare s=%7.3f;" % (self.distance, self.distance)

        ret += "camera {\n"
        ret += "\tlocation <%7.3f, %7.3f, %7.3f>\n" % (1.3 * self.distance, 1.3 * self.distance, -1.3 * self.distance)
        ret += "\tlook_at  <%7.3f, %7.3f, %7.3f>\n" % tuple(0.5 * sum(self.structure.cell[:]))
        ret += "}\n\n"

        if self.structure.nsites > 0:
            d = self.distance
            ret += "light_source { <%7.3f, %7.3f, %7.3f> color White}\n" % (2 * d, 2 * d, 2 * d)

        for imagx in np.arange(-1, 2):
            for imagy in np.arange(-1, 2):
                for imagz in np.arange(-1, 2):

                    for site in self.structure:
                        for symbol in site.symbols:
                            cell = self.structure.cell
                            x = site.position[0] - imagx * cell[0, 0] - imagy * cell[1, 0] - imagz * cell[2, 0]
                            y = site.position[1] - imagx * cell[0, 1] - imagy * cell[1, 1] - imagz * cell[2, 1]
                            z = site.position[2] - imagx * cell[0, 2] - imagy * cell[1, 2] - imagz * cell[2, 2]
                            if (x - self.distance) ** 2 + (y - self.distance) ** 2 + (z + self.distance) ** 2 < 2:
                                continue
                            cr = 0.5 * covalent_radius(symbol)
                            rgb = cpk_colors[atomic_number(symbol)]
                            color = 'rgb < %7.3f, %7.3f, %7.3f>' % (rgb[0], rgb[1], rgb[2])
                            ret += "sphere {\n"
                            ret += "\t<%7.3f, %7.3f, %7.3f>, %7.3f\n\ttexture {\n" % (x, y, z, cr)
                            ret += "\t\tpigment { color %s filter 0.4 transmit %7.3f}\n" % \
                                   (color, 1 - 0.9 * np.exp(-0.1 * (abs(imagx) + abs(imagy) + abs(imagz))))
                            ret += "\t\tnormal { bumps 0.8 scale 0.1 }\n\t\tfinish { phong %7.3f }\n\t}\n}\n\n" % \
                                   np.exp(-0.1 * (abs(imagx) + abs(imagy) + abs(imagz)))

                            if self.structure.nsites <= 0:
                                ret += "light_source { <%7.3f, %7.3f, %7.3f> color White}\n" % (x, y, z)

        ret += """union{
#include "cell.pov"
    scale 1
    rotate <0, 0, 0>
    pigment{rgb <0.3,0.3,0.9>} finish{phong 0.9 ambient 0.42 reflection 0.1}
}
"""
        return ret

    def write_povray(self, filename):

        wf = open(filename, 'w')
        wf.write(self.create_pov())
        wf.close()

        if self.structure.is_crystal:
            self.write_cell('cell.pov')

    def write_cell(self, filename):
        wf = open(filename, 'w')
        ret = ''
        for i in range(3):
            for j in range(3):
                ret += "cylinder { "
                if i == j:
                    ret += " <%7.3f, %7.3f, %7.3f>, " % (0.0, 0.0, 0.0)
                    ret += " <%7.3f, %7.3f, %7.3f>, " % tuple(self.structure.cell[j])
                else:
                    ret += " <%7.3f, %7.3f, %7.3f>, " % tuple(self.structure.cell[i])
                    ret += " <%7.3f, %7.3f, %7.3f>, " % tuple(self.structure.cell[i] + self.structure.cell[j])

                ret += " %7.3f }\n" % (self.distance / 100.0)

            ret += "cylinder { "
            ret += " <%7.3f, %7.3f, %7.3f>, " % tuple(sum(self.structure.cell[:]))
            ret += " <%7.3f, %7.3f, %7.3f>, " % tuple(sum(self.structure.cell[:]) - self.structure.cell[i])

            ret += " %7.3f }\n" % (self.distance / 100.0)
        wf.write(ret)
        wf.close()

# ret += "\topen           // Remove end caps\n"
#                #ret += "\ttexture { %s }\n" % ('T_Stone25 scale 4')
#                ret += "\ttexture { %s }\n" % ('pigment { Col_Glass_Old }')
