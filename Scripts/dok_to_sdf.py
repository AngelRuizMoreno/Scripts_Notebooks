from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS

import sys

import warnings

warnings.filterwarnings("ignore")


def PDBMolSupplier (file=None,sanitize=True):
	mols=[]
	with open(file, 'r') as f:
		doc=[line for line in f.readlines()]
	doc=[line.replace(line.split()[2],line.split()[2].upper()) if 'ATOM' in line else line for line in doc]
	
	scores=[line.split()[-2] for line in doc if 'REMARK Cluster' in line]
	poses=[line.split()[2] for line in doc if 'REMARK Cluster' in line]

	start=[index for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
	finish=[index-1 for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
	finish.append(len(doc))

	interval=list(zip(start,finish[1:]))
	for num,i in enumerate(interval):
		try:
			block = ",".join(doc[i[0]:i[1]]).replace(',','')
			m= Chem.MolFromPDBBlock(block,sanitize=sanitize,removeHs=True)
			m.SetProp('Score',scores[num])
			m.SetProp('Pose',poses[num])
			mols.append(m)
		except Exception:
			pass
	return(mols)

def dok_to_sdf (ref_mol2=None, PDB=None,out_file=None):
	ref=Chem.MolFromMol2File(ref_mol2,removeHs=True)
	mols=PDBMolSupplier(PDB)
	output=Chem.SDWriter(out_file)
	for mol in mols:
		try:
			#newMol = AllChem.AssignBondOrdersFromTemplate(ref, mol)
			#output.write (newMol)

			Chem.Kekulize(ref)
			if mol.GetNumAtoms()==ref.GetNumAtoms():
				atom_pair=[(atom.GetBeginAtomIdx(),atom.GetEndAtomIdx()) for atom in mol.GetBonds()]
				for atom_indx_i, atom_index_j in atom_pair:
					mol.GetBondBetweenAtoms(atom_indx_i,atom_index_j).SetBondType(ref.GetBondBetweenAtoms(atom_indx_i,atom_index_j).GetBondType())
			
			output.write(mol)
		except Exception:
			pass

	output.close()

if __name__ == '__main__':
	print('''
		USAGE: dok_to_sdf.py ref.mol2 docking_results.dok output.sdf

		#NOTE: This script also works to convert .pdb files to aromatic (kekule) .sdf'''
		)
	
	dok_to_sdf(ref_mol2=sys.argv[1], PDB=sys.argv[2],out_file=sys.argv[3])