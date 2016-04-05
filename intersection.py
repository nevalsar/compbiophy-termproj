#CBP Course
#Optimizing short contacts in protein
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
def translate_atom(start_idx, end_idx, atom):

    start = atom_array[start_idx]
    end = atom_array[end_idx-1]


    vector.x_coordinate = end.x_coordinate - start.x_coordinate
    vector.y_coordinate = end.y_coordinate - start.y_coordinate
    vector.z_coordinate = end.z_coordinate - start.z_coordinate

    atom.x_coordinate += vector.x_coordinate
    atom.y_coordinate += vector.y_coordinate
    atom.z_coordinate += vector.z_coordinate

    atom.number -= (end_idx - start_idx)

    return atom

def translate_AtomArray(start_idx, end_idx, atom_array):
    for atom in atom_array:
        if atom.number >= end_idx:
            atom = translate_atom(atom_array[start_idx], atom_array[end_idx-1], atom)

    return atom_array
    
def intersect(m):
	size=len(protein)
	for latom in protein[0:m]:
		for ratom in protein[m:size]:
			if(latom.distance(ratom)<2):
				return true
	return false

def rotate(a,b,c,p):
	size=len(protein)
	pivot=protein[p]
	for ratom in protein[p:size]:
		x=ratom.x-pivot.x
		y=ratom.y-pivot.y
		z=ratom.z-pivot.z
		ratom.x=(cos(a)*cos(b))*x+(cos(a)*sin(b)*sin(c)-sin(a)*cos(c))*y+(cos(a)*sin(b)*cos(c)+sin(a)*sin(c))*z
		ratom.y=(sin(a)*cos(b))*x+(sin(a)*sin(b)*sin(c)+cos(a)*cos(c))*y+(sin(a)*sin(b)*cos(c)-cos(a)*sin(c))*z
		ratom.z=(-sin(b))*x+(cos(b)*sin(c))*y+(cos(b)*cos(c))*z
		ratom.x=ratom.x+pivot.x
		ratom.y=ratom.y+pivot.y
		ratom.z=ratom.z+pivot.z

def findangles():
	a=0#Roll
	b=0#Pitch
	c=0#Yaw
	readprotein(atoms)
	translate(atoms,st,en)
	batoms=atoms
	while a<180:
		while b<180:
			while c<180:
				rotate(a,b,c)
				if(intersect(m)==false):
					print 'Found Valid Rotation for no short contacts:'
					print 'Roll,Pitch,Yaw:'
					print a,b,c

			c=c+delta
		b=b+delta
	a=a+delta
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
scut=10
ecut=20
for patom in protein:
	for ratom in rna:
		if patom.distance(ratom)<2:
			try:
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)][patom.label][ratom.label]+=1
			except:
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)]={'Sidechain':{'Sugar':0,'Phosphate':0,'Base':0},'Backbone':{'Sugar':0,'Phosphate':0,'Base':0}}
				interaction[patom.chain_id+':'+patom.residue_name+'_'+str(patom.residue_no)][patom.label][ratom.label]=1