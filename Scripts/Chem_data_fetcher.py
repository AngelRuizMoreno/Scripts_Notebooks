import requests
import argparse
from bs4 import BeautifulSoup
import sys

def lookup (identifier=None,out=None):
	try:
		url = 'https://cactus.nci.nih.gov/chemical/structure/' + identifier + '/' +out
		with requests.get(url) as response:
			http = response.content
			result = BeautifulSoup(http,"html.parser")
		return result
	except Exception:
		pass

def main():
	print('''
Outputs available:

Satandard InChIKey       = stdinchikey
Standard InChI           = stdinchi
SMILES                   = smiles
Molecular Weight         = mw
Violations to Rule of 5  = rule_of_5_violation_count
Names                    = names
CAS number               = cas
Chemical formula         = formula
IUPAC name               = iupac_name
''')
	
	input_file=sys.argv[1]
	output_file=sys.argv[2]
	representation=sys.argv[3]

	out=open(output_file,'w', encoding="utf-8")
	with open(input_file,'r',encoding="utf-8") as f:
		for line in f:
			try:
				if 'FIELD CA Index Name:' in line:
					line=line.replace('FIELD CA Index Name:','').split(', ')[0:2]
					line=''.join([i.rstrip(',').replace('\n','') for i in reversed(line)])
					r=lookup(line,representation)
					out.write(r.string+' {} \n'.format(line))
			except Exception:
				pass
	out.close()

if __name__ == '__main__':
	main()