from rdkit import Chem

def find_connected_Hs(atom_idx,mol_w_H):
    '''
    :param atom_idx: index of atom in question
    :param mol_w_H: apply Chem.AddHs(mol) to it
    :return: Indices of all H-atoms connected to an Atom
    '''
    H_indices = []
    atom = mol_w_H.GetAtomWithIdx(atom_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:  # Wasserstoffatom
            H_indices.append(neighbor.GetIdx())
    return H_indices

def find_substruc(mol: Chem.Mol, submol):
    '''
    Finds the Indices of all atoms (includings Hs) that belong to a substructure in a molecule
    :param mol:
    :param submol: mol (from smart or not) of substruc contained in molecule
    :return: indices
    '''
    mol_w_H = Chem.AddHs(mol)
    match = mol_w_H.GetSubstructMatch(submol)
    H_idx = []
    for atom in list(match):
        H_idx += find_connected_Hs(atom,mol_w_H)

    match = list(match) + H_idx
    return match

def mol_mapper(mol: Chem.Mol,substruc: Chem.Mol,linker: Chem.Mol):
    mol_H = Chem.AddHs(mol)
    Tail = find_substruc(mol, substruc)
    Leaving_group = list(set(range(mol_H.GetNumAtoms())) - set(Tail))
    Linker = find_substruc(mol_H,linker)
    return Tail, Leaving_group, Linker

def linker_atom_mapper(linker_idx,mol: Chem.Mol,Tail_idx,LG_idx):
    '''
    :param linker_idx:
    :param mol:
    :param Tail_idx:
    :param LG_idx:
    :return: Returns indexes of several interesting atom positions
    '''
    mol_H = Chem.AddHs(mol)
    for idx in linker_idx:
        # find first N and attached H
        if idx in Tail_idx and mol_H.GetAtomWithIdx(idx).GetAtomicNum() == 7:
            N_1 = idx
            for neighbour in mol_H.GetAtomWithIdx(N_1).GetNeighbors():
                if neighbour.GetAtomicNum() == 1:
                    H_at_N = neighbour.GetIdx()

        #get carbonyl O and C
        elif idx in Tail_idx and mol_H.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            carbonyl_O = idx

        elif idx in Tail_idx and mol_H.GetAtomWithIdx(idx).GetAtomicNum() == 6:
            carbonyl_C = idx

        #get linking atom of LG
        elif idx in LG_idx and mol_H.GetAtomWithIdx(idx).GetAtomicNum() in [7,8]:
            linker_atom_idx = idx

            #get neighbors of linking atom of LG that are part of LG
            Neighbours_linker_LG = []
            for neighbour in mol_H.GetAtomWithIdx(linker_atom_idx).GetNeighbors():
                if neighbour.GetIdx() in LG_idx:
                    Neighbours_linker_LG.append(neighbour.GetIdx())

    return N_1, H_at_N, carbonyl_O, carbonyl_C, linker_atom_idx, Neighbours_linker_LG





if __name__ == "__main__":
    smile_phenol = 'OC([C@@H](NC(C)=O)CCCNC(OC1=CC=CC=C1)=O)=O'
    smile_imidazol = 'OC([C@@H](NC(C)=O)CCCNC(N1C=CN=C1)=O)=O'
    smile_triazol = 'OC([C@@H](NC(C)=O)CCCNC(N1C=NN=C1)=O)=O'
    smile_pyridin = 'OC([C@@H](NC(C)=O)CCCNC(NC1=C(C)C=CC=N1)=O)=O'

    # substruc = 'OC([C@@H](NC(C)=O)CCCNC(*)=O)=O'
    substruc = 'OC([C@@H](NC(C)=O)CCCNC=O)=O'
    linker_smart = '[N&H1](C([O,N,n])=O)'

    mol_smi = Chem.MolFromSmiles(smile_pyridin)

    T, LG, Li = mol_mapper(mol_smi,Chem.MolFromSmarts(substruc), Chem.MolFromSmarts(linker_smart))

    N, H, O, C, L, nei = linker_atom_mapper(Li, mol_smi,T,LG)

    from Utils import draw_indexed_molecule
    draw_indexed_molecule(mol_smi,'/home/student/j_spie17/molecular_prosthetics/test/','Pyridin')




    print('done')