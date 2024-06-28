from rdkit import Chem
from morfeus import *
import numpy as np

#this is incomplete!!!! only important elements in here yet
periodic_table = {1 : 'H', 2 : 'He' , 3 : 'Li', 4 : 'Be', 5 : 'B', 6 : 'C' , 7 : 'N', 8 : 'O' , 9 : 'F', 10 : 'Ne', 11 : 'Na',
                      12 : 'Mg', 13 : 'Al', 14 : 'Si', 15 : 'P', 16 : 'S', 17 : 'Cl', 35 : 'Br', 53 : 'I'}


class Structure(object):
    '''
    Stores the geometrical Structure of a Smile from coordinates and relates it to a rdkit mol object
    '''

    def __init__(self, name=None, smile=None, elements=None, coordinates=None, geom_indexes=None):
        self.name = name
        self.smile = smile

        # Geometry related
        self.geom_path = None
        self.elements = elements
        self.coordinates = coordinates
        self.geom_indexes = geom_indexes

        #rdkit related
        self.mol = None


    def read_gauss_output(self, log_path):
        '''
        Read a Gaussian log file and extract the coordinates of the optimized Geometry, as well as Elements, Indexes and Energy
        :param log_path: .log file
        :return: Energy, Indexes, Elements, Coordinates of each atom
        '''

        #read in file
        self.geom_path = log_path
        with open(self.geom_path, 'r') as f:
            lines = f.readlines()
        lines_cleaned = [line.strip() for line in lines]

        #Extract Geometry
        '''
        parser that iterates though each line of the script, when a 'standard orientation' is found the following coordinates are saved, 
        all sets of coordinates are saved in standard_orientation_block
        '''
        standard_orientation_block = []
        standard_orientation = []
        inside_block = False
        count_hyphens = 0


        for line in lines_cleaned:
            if "Standard orientation" in line:
                inside_block = True
            elif inside_block:
                if "-----" in line:
                    count_hyphens += 1
                    if count_hyphens == 3:
                        inside_block = False
                        standard_orientation_block.append(standard_orientation)
                        standard_orientation = []
                        count_hyphens = 0
                else:
                    if count_hyphens == 2:
                        standard_orientation.append(line)

        #Convert text output into a coordinates matrix and a list of elements
        orientation = standard_orientation_block[-1]  # last geometry(orientatiion) should be the optimized one
        orientation = [line.split() for line in orientation]  # remove spaces
        array = np.array(orientation)
        self.coordinates = array[:, 3:].astype(float)
        elements_idx = array[:, 1].astype(int)

        self.elements = [periodic_table[key] for key in elements_idx] #Convert atomic numbers to element symbol ( 1 -> 'H')

        return




if __name__ == "__main__":
    smile = "C(C)CSC(=O)c1ccccc1"
    Strx = Structure(smile=smile)
    Strx.name = 'struc1'




    print('done')