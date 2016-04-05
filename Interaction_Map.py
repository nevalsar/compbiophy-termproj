#CSB Course
#Assignment 2
#Protein RNA Complex Interaction of 1ASY.pdb
#Shrihari Bhat
#11CS30008

class Atom:
	def __init__(self,line):
		#Sample line of ATOM:
		#0         1         2         3         4         5         6         7       
		#01234567890123456789012345678901234567890123456789012345678901234567890123456789 
		#ATOM   6338  CA  TYR A 456      60.059  70.952 -26.972  1.00 23.16           C
		self.srno=int(line[6:11])
		self.name=line[12:16].strip()
		self.residue_name=line[17:20].strip()
		self.residue_no=int(line[22:26])
		self.chain_id=line[21]
		self.x=float(line[30:38])
		self.y=float(line[38:46])
		self.z=float(line[46:54])
		self.data=line
		#Decide the Label of atom from RNA(Sugar/Phosphate/Base) OR Protein(Backbone(N/CA/C)/Sidechain)
		if len(self.residue_name)<3:
			if "'" in self.name:
				self.label='Sugar'
			elif "P" in self.name:
				self.label='Phosphate'
			else:
				self.label='Base'
		else:
			if self.name in ['N','CA','C']:
				self.label='Backbone'
			else:
				self.label='Sidechain'

	def distance(self,other):
		xs=self.x-other.x
		ys=self.y-other.y
		zs=self.z-other.z
		square_dist=xs*xs+ys*ys+zs*zs
		dist=square_dist**0.5
		return dist

protein=[]
rna=[]
# dic={}
pdb=open('1ASY.pdb','r')
# print 'File Opened'
for line in pdb:
	if line.startswith('ATOM'):
		a=Atom(line)
		# try:
		# 	dic[line[21]]+=1
		# except:
		# 	dic[line[21]]=0
		# 	dic[line[21]]=1
		if a.chain_id=='B':
			# print 'Protein'
			#Protein Atom
			protein.append(a)
		if a.chain_id=='S':
			#RNA Atom
			# print 'RNA'
			rna.append(a)

interaction={}
for patom in protein:
	for ratom in rna:
		if patom.distance(ratom)<4.5:
			try:
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)][patom.label][ratom.label]+=1
			except:
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)]={'Sidechain':{'Sugar':0,'Phosphate':0,'Base':0},'Backbone':{'Sugar':0,'Phosphate':0,'Base':0}}
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)][patom.label][ratom.label]=1

print '{0:>28s}{1:>13s}{2:>9s}'.format('Sugar','Phosphate','Base')
for residue in interaction:
	for key1 in ['Sidechain','Backbone']:
		print '{0:<10s}{1:>10s}'.format(residue,key1),'{0:^10d}'.format(interaction[residue][key1]['Sugar']),'{0:^10d}'.format(interaction[residue][key1]['Phosphate']),'{0:^10d}'.format(interaction[residue][key1]['Base'])
	print