from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import argparse

def ConformersFromSmiles(smiles, nConfs=10, maxIters=500, verbose=False):
    '''
    Generate conformers for a molecule specified by its SMILES string.

    Parameters:
    - smiles (str): SMILES representation of the molecule.
    - nConfs (int, optional): Number of conformers to generate. Default is 10.
    - maxIters (int, optional): Maximum number of force-field optimization iterations for each conformer. Default is 500.
    - verbose (bool, optional): If True, print information about conformer minimization convergence. Default is False.

    Returns:
    Chem.Mol: The molecule with generated conformers.
    '''
    mol = Chem.MolFromSmiles(smiles)
    mol = AllChem.AddHs(mol)

    confs = AllChem.EmbedMultipleConfs(mol, nConfs, randomSeed=42, 
                                        useRandomCoords=False,
                                        enforceChirality=True,
                                        useExpTorsionAnglePrefs=True, 
                                        useBasicKnowledge=True,
                                        useSmallRingTorsions=True, 
                                        useMacrocycleTorsions=True,  
                                        )
    converged = 0
    for conf in confs:
        opt = AllChem.MMFFOptimizeMolecule(mol, confId=conf, maxIters=maxIters)
        converged+=opt
        if (opt == -1) and verbose:
            print(f'Force field could not be setup for conformer {conf}!')
        else:
            converged += opt
    if verbose:
        print(f'{converged} conformer minimisations failed to converge')
    return mol

def ConfsToAlignedMolsList(multiConfMol):
    '''
    Align conformers within a molecule and convert them to a list of aligned molecules.

    Parameters:
    - multiConfMol (Chem.Mol): A molecule with multiple conformers.

    Returns:
    List[Chem.Mol]: A list of aligned molecules derived from the conformers.
    '''
    AllChem.AlignMolConformers(multiConfMol)
    mols = []
    cids = [x.GetId() for x in multiConfMol.GetConformers()]
    for cid in cids:
        mol = Chem.MolToMolBlock(multiConfMol, confId=cid)
        mol = Chem.MolFromMolBlock(mol, removeHs=False)
        mols.append(mol)
    return mols

def ComputeMedianConformer(mol):
    '''
    Compute the median conformer of a molecule as the conformer whose atomic coordinates 
    are closest to the mean coordinates of all conformers.
    
    Parameters:
    - mol (Chem.Mol): A molecule with one or more conformers.
    
    Returns:
    Chem.Mol: A new molecule with a single conformer, which is the median conformer.
    '''
    confs = mol.GetConformers()
    if not confs:
        raise ValueError('No conformers found in the molecule.')
    
    num_atoms = mol.GetNumAtoms()
    nConfs = len(confs)
    
    mean_coords = []
    for i in range(num_atoms):
        x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
        for conf in confs:
            pos = conf.GetAtomPosition(i)
            x_sum += pos.x
            y_sum += pos.y
            z_sum += pos.z
        mean_x = x_sum / nConfs
        mean_y = y_sum / nConfs
        mean_z = z_sum / nConfs
        mean_coords.append((mean_x, mean_y, mean_z))
    
    best_conf = None
    best_dist = float('inf')
    for conf in confs:
        total_sq_dist = 0.0
        for i in range(num_atoms):
            pos = conf.GetAtomPosition(i)
            dx = pos.x - mean_coords[i][0]
            dy = pos.y - mean_coords[i][1]
            dz = pos.z - mean_coords[i][2]
            total_sq_dist += dx*dx + dy*dy + dz*dz
        if total_sq_dist < best_dist:
            best_dist = total_sq_dist
            best_conf = conf
    
    if best_conf is None:
        raise RuntimeError('Could not determine a median conformer.')

    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_conf = Chem.Conformer(num_atoms)

    for i in range(num_atoms):
        pos = best_conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(i, pos)
    new_mol.AddConformer(new_conf, assignId=True)
    
    return new_mol

def main():
    parser = argparse.ArgumentParser(description='Generate the mean conformer from a SMILES string and output as a PDB file.')
    parser.add_argument('-i', '--input', required=True, help='Input SMILES string')
    parser.add_argument('-o', '--output', required=True, help='Output file name (without .pdb extension)')
    parser.add_argument('--nconfs', type=int, default=10, help='Number of conformers to generate (default: 10)')
    parser.add_argument('--maxiters', type=int, default=500, help='Maximum force-field optimization iterations per conformer (default: 500)')
    parser.add_argument('--verbose', action='store_true', help='Print detailed output during processing')
    
    args = parser.parse_args()
    
    mol = ConformersFromSmiles(args.input, nConfs=args.nconfs, maxIters=args.maxiters, verbose=args.verbose)
    mean_mol = ComputeMedianConformer(mol)
    pdb_filename = 'outputs/' + args.output + '.pdb'
    Chem.MolToPDBFile(mean_mol, pdb_filename)
    print(f'Mean conformer saved to {pdb_filename}')

if __name__ == '__main__':
    main()