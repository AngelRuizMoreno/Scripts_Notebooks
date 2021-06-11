import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from itertools import groupby
import math
import sys


def analyzer(file=None):
	print ('Loading results')
	table=pd.read_excel(file,index_col=[0],engine='openpyxl')

	interaction_types=list(set(table['Type']))

	new_names=[chain+':'+resname for chain,resname in zip(table['Chain'],table['Residue'])]
	table['newNames']=new_names
	heatmap=pd.DataFrame(index=sorted(list(set(table['newNames']))), columns=sorted(interaction_types))

	print ('initializing, saving data and ploting results')

	for i in list(set(table['newNames'])):
		residue=sorted(list(table[table['newNames']==i]['Type']))
		groups=groupby(residue)
		for x in groups:
			heatmap.loc[i,x[0]]=len(list(x[1]))
	
	heatmap.to_excel('Heatmap_summary.xlsx')

	plot_heatmap(data=heatmap)
	get_frequency(table=table)

def plot_heatmap(data=None):

	fig, ax = plt.subplots(figsize=(10,2))
	
	ax=sns.heatmap(data.transpose().fillna(0),cmap='Blues',cbar_kws=dict(label='Frequency',shrink=1,orientation='vertical',spacing='uniform',pad=0.02))

	plt.title('Interactions Summary',size='22',weight='bold')
	plt.xlabel('Residue',fontsize=20,fontweight='bold')
	plt.xticks (rotation=90,fontsize=12)
	plt.yticks (fontsize=12)
	plt.tick_params ('both',width=2,labelsize=12)
	plt.savefig('HeatMap.png',dpi=300,format='png', bbox_inches = 'tight')

def get_frequency(table=None):
	time=table['Time'].drop_duplicates()
	index=table.index.drop_duplicates()
	data=pd.DataFrame()
	
	for x in index: 
		data.loc[x,'H-bon']=len([i for i in table.loc[x,'Type'] if i=='H-bond'])
		data.loc[x,'Water-bridge']=len([i for i in table.loc[x,'Type'] if i=='Water-bridge'])
		data.loc[x,'Salt-bridge']=len([i for i in table.loc[x,'Type'] if i=='Salt-bridge'])
		data.loc[x,'Pi-cation']=len([i for i in table.loc[x,'Type'] if i=='Pi-cation'])
		data.loc[x,'Hydrophobic']=len([i for i in table.loc[x,'Type'] if i=='Hydrophobic'])
		data.loc[x,'Pi-stacking']=len([i for i in table.loc[x,'Type'] if i=='Pi-stacking'])
		data.loc[x,'X-bond']=len([i for i in table.loc[x,'Type'] if i=='X-bond'])
		data.loc[x,'Metal-complex']=len([i for i in table.loc[x,'Type'] if i=='Metal-complex'])
	data.to_excel('Frequency.xlsx')
	
	return plot_frequency(data=data)

def plot_frequency(data=None):
	
	for column in data.columns:
		
		plt.rcParams['axes.linewidth'] = 1.5
		fig, ax = plt.subplots(figsize=(10,5))

		y_values=[float('nan') if x==0 else x for x in data[column]]

		ax.scatter(data.index,y_values,c='navy',s=1, marker='.')

		plt.title (column,fontsize=24,fontweight='bold')
		plt.xlabel ('Frame',fontsize=20,fontweight='bold')
		plt.ylabel ('Frequency',fontsize=20,fontweight='bold')

		yint = range(math.floor(min(data[column])), math.ceil(max(data[column]))+1)

		plt.yticks(yint)

		plt.tick_params ('both',width=2,labelsize=14)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		plt.tight_layout()
		plt.savefig(column+'.png',dpi=300,format='png',bbox_inches='tight')

if __name__ == "__main__":
	print ('USAGE: python plipMD_analizer.py Results.xlsx')
	analyzer(sys.argv[1])