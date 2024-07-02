from rdkit import Chem
from rdkit.Chem import Draw
from morfeus import *
import numpy as np
import os
from Utils import  *
from rdkit.Chem import rdDetermineBonds

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

            #now check for thermal Free Energy
            if inside_block == False and "Sum of electronic and thermal Free Energies=" in line:
                Energy = float(line.split()[-1])
                self.energy = Energy * 627.15 #Energy in kcal/mol
            else:
                continue

        #Convert text output into a coordinates matrix and a list of elements
        orientation = standard_orientation_block[-1]  # last geometry(orientatiion) should be the optimized one
        orientation = [line.split() for line in orientation]  # remove spaces
        array = np.array(orientation)
        self.coordinates = array[:, 3:].astype(float)
        elements_idx = array[:, 1].astype(int)

        self.elements = [periodic_table[key] for key in elements_idx] #Convert atomic numbers to element symbol ( 1 -> 'H')

        return

    def rdkit_geometries(self):
        self.mol = Chem.MolFromSmiles(self.smile)
        return

    def draw_indexed_molecule(self, save_dir=None):
        '''
        Saves an image of the Structure with H_atoms and all Indexes
        :param save_dir: Location where the .png will be saved
        :return:
        TODO: As of yet the directory is not changed back from the save_dir
        '''
        if not save_dir:
            save_dir = input('Please set a path to save the .png file:')
        os.chdir(save_dir)

        if not self.name:
            self.name = input('Please set a name to save the .png file:')

        if not self.mol:
            self.mol = Chem.MolFromSmiles(self.smile)

        self.mol_H = Chem.AddHs(self.mol)

        # Nummeriere alle Atome einschließlich Wasserstoff
        for atom in self.mol_H.GetAtoms():
            atom.SetProp('atomLabel', str(atom.GetIdx()))
        # Zeichne die Molekülstruktur und speichere das Bild

        img = Draw.MolToImage(self.mol_H, size=(300, 300), kekulize=True, wedgeBonds=True)
        img.save(f"{self.name}_indexed.png")


    def coordinates_to_mol_via_xyz(self,tmp_dir,delete_xyz=True):
        '''
        This function (temporarily) creates an xyz file from previous read in coordinates and elements
        This in turn is read in with rdkit function to create a Mol object that has the same Indexing as the xyz file and has the geometries in its conformer
        :param tmp_dir: Where the .xyz file will be created
        :param delete_xyz: optinal if xyz file should be deleted afterwards
        :return:
        '''
        #change directory
        if not tmp_dir:
            tmp_dir = input('Please set a path to save the .xyz file:')
        os.chdir(tmp_dir)

        if not self.name:
            self.name = input('Please set a name to save the .xyz file:')

        #save geometry as xyz
        write_xyz_file(filename=f"{self.name}.xyz",atoms=self.elements,coordinates=self.coordinates)

        #read xyz
        raw_mol = Chem.MolFromXYZFile(f"{self.name}.xyz")
        conn_mol = Chem.Mol(raw_mol)

        rdDetermineBonds.DetermineConnectivity(conn_mol, charge=0)
        rdDetermineBonds.DetermineBondOrders(conn_mol, charge=0)
        # Chem.RemoveHs(conn_mol) #TODO: should they stay or should they go?
        self.mol = conn_mol

        #update smiles
        self.smile = Chem.CanonSmiles(Chem.MolToSmiles(self.mol))

        #delete xyz
        if delete_xyz == True:
            os.remove(f"{self.name}.xyz")

        return
    




if __name__ == "__main__":
    smile = "C(C)CSC(=O)c1ccccc1"
    Strx = Structure(smile=smile)
    Strx.name = 'struc1'

    Strx.read_gauss_output("/home/student/j_spie17/molecular_prosthetics/gaussian/YPACPXJQEYROHD-UHFFFAOYSA-N_conf_0.log")


    from rdkit.Chem import rdDetermineBonds

    raw_mol = Chem.MolFromXYZFile("/home/student/j_spie17/molecular_prosthetics/test/YPACPXJQEYROHD-UHFFFAOYSA-N_conf_0.xyz")
    conn_mol = Chem.Mol(raw_mol)

    rdDetermineBonds.DetermineConnectivity(conn_mol, charge=0)
    rdDetermineBonds.DetermineBondOrders(conn_mol, charge=0)
    # Chem.RemoveHs(conn_mol)

    print(Chem.MolToSmiles(conn_mol))



    print('done')