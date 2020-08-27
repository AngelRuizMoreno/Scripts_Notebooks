import sys
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, Draw,rdMolDescriptors

def get_descriptors (mols=None):
	table=pd.DataFrame()
	i=0
	print ('Computing descriptors for {} molecules'.format(len(mols)))
	for mol in mols:
		if mol:
			table.loc[i,'Smiles']=Chem.MolToSmiles(mol,True)
			table.loc[i,'MolWt']=Descriptors.MolWt(mol)
			table.loc[i,'LogP']=Descriptors.MolLogP(mol)
			table.loc[i,'NumHAcceptors']=Descriptors.NumHAcceptors(mol)
			table.loc[i,'NumHDonors']=Descriptors.NumHDonors(mol)
			table.loc[i,'NumHeteroatoms']=Descriptors.NumHeteroatoms(mol)
			table.loc[i,'NumRotableBonds']=Descriptors.NumRotatableBonds(mol)
			try:
				AllChem.EmbedMolecule(mol,useRandomCoords=True)
				AllChem.MMFFOptimizeMolecule(mol,mmffVariant='MMFF94s')
				table.loc[i,'TPSA']=Descriptors.TPSA(mol)
				table.loc[i,'PMI1']=Descriptors3D.PMI1(mol)
				table.loc[i,'PMI2']=Descriptors3D.PMI2(mol)
				table.loc[i,'PMI3']=Descriptors3D.PMI3(mol)
				table.loc[i,'PBF']=rdMolDescriptors.CalcPBF(mol)
				table.loc[i,'NPR1']=rdMolDescriptors.CalcNPR1(mol)
				table.loc[i,'NPR2']=rdMolDescriptors.CalcNPR2(mol)
				table.loc[i,'ISF']=Descriptors3D.InertialShapeFactor(mol)
				i=i+1
			except Exception:
				i=i+1
				continue

	table.to_csv ('Descriptors.csv')
	print ('Computed descriptors in file Descriptors.csv')           
	return (plot_NPR(table=table))
	

def plot_NPR(table=None):

	plt.rcParams['axes.linewidth'] = 1.5
	plt.figure(figsize=(10,10))

	ax=sns.scatterplot('NPR1','NPR2',data=table,s=50,palette=sns.color_palette("husl", 1),linewidth=0.25)
	x1, y1 = [0.5, 0], [0.5, 1]
	x2, y2 = [0.5, 1], [0.5, 1]
	x3, y3 = [0,1],[1,1]
	plt.plot(x1, y1,x2,y2,x3,y3,c='gray',ls='--',lw=2)

	plt.xlabel ('NPR1',fontsize=12,fontweight='bold')

	plt.ylabel ('NPR2',fontsize=12,fontweight='bold')

	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

	#plt.legend(loc='lower left',frameon=False,prop={'size': 12},ncol=1)
	plt.tick_params ('both',width=2,labelsize=10)
	plt.tight_layout()
	print ('Rendering Principal Moment of Inertia plot')
	plt.savefig('NPR.png',dpi=300,quality=95,format='png')
	print ('All done, have fun')
	#plt.show()

def main():
		if '.sdf' in sys.argv[1]:
			mols=Chem.SDMolSupplier(sys.argv[1])
			get_descriptors(mols=mols)
		
		if '.smiles' in sys.argv[1]:
			mols=Chem.SmilesMolSupplier(sys.argv[1])
			get_descriptors(mols=mols)

if __name__ == "__main__":
	main()