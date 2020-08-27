from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem, rdMolAlign
import sys


def get_conformers(smiles=None, anchor=None, num_confs=None, output=None):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(mol)

    constrain = Chem.SDMolSupplier(anchor,False)[0]

    r = rdFMCS.FindMCS([mol, constrain])
    a = mol.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    b = constrain.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    amap = list(zip(a, b))

    coors = dict()

    for i in a:
        coors[i] = mol.GetConformer().GetAtomPosition(i)

    w = Chem.SDWriter(output)

    AllChem.EmbedMolecule(mol)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    for i in mol.GetSubstructMatch(constrain):
        ff.MMFFAddPositionConstraint(i, 0, 1.0e5)

    confs = AllChem.EmbedMultipleConfs(mol,
                                       numConfs=int(num_confs),
                                       pruneRmsThresh=0.75,
                                       coordMap=coors,
                                       enforceChirality=True,
                                       useExpTorsionAnglePrefs=True,
                                       useBasicKnowledge=True)

    for element in confs:
        rmsd = rdMolAlign.AlignMol(mol, constrain, element, 0, atomMap=amap)
        w.write(mol, confId=element)
    w.close()


def main():
    print('''
        USAGE: python Constrained_confs.py smiles_code anchor.sdf num_confs output.sdf
        ''')

    get_conformers(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == '__main__':
    main()
