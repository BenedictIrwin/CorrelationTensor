from numpy import zeros, save, sqrt, vectorize, array, transpose, corrcoef, mean, dot, matmul, subtract, load, correlate, vectorize
from numpy.linalg import eigh, inv
from os.path import isfile
from scipy.stats.stats import pearsonr

NN = 12852 

## Load all of the data in
Mus = load("GIANTModel_Mus.npy")
#Atoms = load("Model_Atoms.npy")
#Strings = load("Model_Strings.npy")
Us = load("GIANTModel_Us.npy")
Trajs = load("GIANTModel_Trajs.npy")

N_atoms = len(Mus)
print(N_atoms)
## Atom a

Ks = zeros(shape=(N_atoms,N_atoms,3,3))
  
for a in range(N_atoms):
  ### Atom b
  for b in range(N_atoms):
    if(b==0): print( "{0: 0.3f}%".format(100*(a*N_atoms + b)/(N_atoms*N_atoms)) )
    correlation_matrix = array([ [ pearsonr(Trajs[a][k],Trajs[b][l])[0] for l in [0,1,2]] for k in [0,1,2] ])
    Ks[a,b] = matmul(matmul(Us[a], correlation_matrix),transpose(Us[b]))

Ks = array(Ks)
print(Ks.shape)
save("GIANTModel_Coeffs",Ks)
exit()


## Define an ordering to pairs of atoms


## For each pair of atoms generate a coefficient







exit()





## To convert a standard pdb string tyo a vector of x y z coords
def pdb_to_vec(s): return [float(s[30:38]),float(s[38:46]),float(s[46:54])]
def to_vec(s): return [pdb_to_vec(i) for i in s]
#pdb_vec = vectorize(pdb_to_vec)

## Take two amino classes and generate an interaction dictionary
## This is a set of coefficients between the atoms
def generate_coefficient_dictionary(a1,a2):
  local_dict = {}

  for i in range(a1.numatoms):
    for j in range(a2.numatoms):
      correlation_matrix = [ [ pearsonr(a1.projections[k][i],a2.projections[l][j])[0] for l in [1,2,3]] for k in [1,2,3] ] 
      key = a1.string+":"+a1.atomnames[i]+":"+a2.string+":"+a2.atomnames[j]
      #key=key.replace(" ","_")
     
      ## U is backwartds! the order is ascending!
      print("U is backwards!!!!")
      exit()
 
      ## The key is going to be 
      key = (a1.atomnums[i]-1) * NN*(a2.atomnums[j]-1)

      ## For these coeffs we should multiply the raw xyz force on the left hand side
      local_dict[key] = matmul(matmul(a1.Us[i], correlation_matrix),a2.Uinvs[j])
  print(local_dict)
  return local_dict

## Make an amino class to be used for all analysis
class Amino:
  def __init__(self,string):
    self.string = string 
    self.start = None
    self.stop  = None
    self.store = []

    ## Store the centroids of each atom in the amino
    self.mus = []
    self.Us = []
    self.Uinvs = []

    ## Store the eigenvalues and vectors of each atom in the vector
    #self.w1s = []
    self.v1s = []
    #self.w2s = []
    self.v2s = []
    #self.w3s = []
    self.v3s = []

    ## Store a projection of motion across all frames to derive correlations
    self.projections = {1 : [], 2 : [], 3 : []}
    self.projections_mat = []

    ## Naming information
    self.atomnames = []
    self.numatoms = None
    
  ## Find the line number in a test file
  def locate(self,filename):
    indices = []
    lines = []
    with open(filename,"r") as f:
      for i, line in enumerate(f):
        if(self.string in line):
          indices.append(i)
          lines.append(line)
    if(len(indices)==0):
      print("Could not find amino {}".format(self.string))
      exit()
    self.atomnames = [ q[13:15] for q in lines  ]
    self.atomnums = [ int(q[6:12]) for q in lines  ]
    self.numatoms = len(self.atomnames)
    self.start = min(indices)
    self.stop = max(indices) + 1 ## Add one for ranges

  def add_data(self,flines):
    ## We read the data from the file
    self.store.append(to_vec(flines[self.start:self.stop]))

  ## Calculate the principle direction of motion
  ## Calculate the size of the eigenvalue
  ## Calculate the projections
  def calculate_essentials(self,f1,f2,f3,f4,f5):

    self.store = transpose( array(self.store), (1,2,0) )
    self.num_atoms = self.store.shape[0]
    for i in self.store:
      self.mus.append(mean(i, axis = 1))
      ## Hermitian matrix
      w, v = eigh(corrcoef(i))
      ## The eigenvalue matrix used to project
      self.Us.append(v[:,::-1])
      self.Uinvs.append(inv(v[:,::-1]))

    ## Calculate projections
    ## Each atom has three projections
    ## Two atoms have 3x3 combinations of projections
    for i in range(self.num_atoms):
      d = transpose(self.store[i])
      mu = self.mus[i]
      qqq = matmul(self.Uinvs[i],transpose(subtract(d,mu)))
      name = self.atomnames[i]

      f1.write("{}\n".format(self.string))
      f2.write("{}\n".format(name))
      f3.write("{}\n".format(mu.tobytes()))
      f4.write("{}\n".format(self.Us[i].tobytes()))
      #f.write("{},".format(self.Uinvs[i].tobytes()))
      f5.write("{}\n".format(qqq.tobytes()))

 
## Define a complete file that represents the structure of the protein for all frames
test_file = "FILE1/Stage1/s0"
equilibrated_point = 25000
end_frame = 25100

## From the structure (manually) define some key residues in the structure
binding_site_amino_strings = ["PHE I1504","ARG I1506","PHE E 773","ARG E 775","TYR A 149","TYR A 197","TYR G1241","TYR G1289"]
top_annulus_amino_strings = ["LYS A 266","LYS C 625","LYS E 987","LYS G1358","LYS I1718"]
middle_annulus_amino_strings = ["ARG A 261","ARG C 620","ARG E 982","ARG G1353","ARG I1713"]
bottom_annulus_amino_strings = ["ARG A 242","ARG C 601","ARG E 963","ARG G1334","ARG I1694"]
OH_annulus_1_strings = ["THR A 263","SER C 622","SER E 984","THR G1355","SER I1715"]
OH_annulus_2_strings = ["THR A 255","THR C 614","THR E 976","THR G1347","THR I1707"]
OH_annulus_3_strings = ["THR A 252","THR C 611","THR E 973","THR G1344","THR I1704"]
OH_annulus_4_strings = ["THR A 249","THR C 608","THR E 970","THR G1341","THR I1701"]
hydrophobic_annulus_1_strings = ["LEU A 264","LEU C 623","LEU E 985","LEU G1356","LEU I1716"] 
hydrophobic_annulus_2_strings = ["ILE A 256","LEU C 615","LEU E 977","ILE G1348","LEU I1708"] 
hydrophobic_annulus_3_strings = ["LEU A 251","LEU C 610","LEU E 972","LEU G1343","LEU I1703"] 
hydrophobic_annulus_4_strings = ["ILE A 247","ILE C 606","VAL E 968","ILE G1339","VAL I1699"] 

## All methionines very interesting!
## This is roughly the place to try and squish down
center_of_channel = ["MET A 253", "MET C 612", "MET E 974", "MET G1345", "MET I1705"]

## Define the aminos as individual objects
binding_site_aminos = [Amino(s) for s in binding_site_amino_strings]
top_annulus_aminos = [Amino(s) for s in top_annulus_amino_strings]
middle_annulus_aminos = [Amino(s) for s in middle_annulus_amino_strings]
bottom_annulus_aminos = [Amino(s) for s in bottom_annulus_amino_strings]
OH_annulus_1_aminos = [Amino(s) for s in OH_annulus_1_strings]
OH_annulus_2_aminos = [Amino(s) for s in OH_annulus_2_strings]
OH_annulus_3_aminos = [Amino(s) for s in OH_annulus_3_strings]
OH_annulus_4_aminos = [Amino(s) for s in OH_annulus_4_strings]
hydrophobic_annulus_1_aminos = [Amino(s) for s in hydrophobic_annulus_1_strings]
hydrophobic_annulus_2_aminos = [Amino(s) for s in hydrophobic_annulus_2_strings]
hydrophobic_annulus_3_aminos = [Amino(s) for s in hydrophobic_annulus_3_strings]
hydrophobic_annulus_4_aminos = [Amino(s) for s in hydrophobic_annulus_4_strings]

## Define a list as a global collection
collection = [binding_site_aminos, top_annulus_aminos, middle_annulus_aminos, bottom_annulus_aminos, OH_annulus_1_aminos, OH_annulus_2_aminos, OH_annulus_3_aminos, OH_annulus_4_aminos, hydrophobic_annulus_1_aminos, hydrophobic_annulus_2_aminos, hydrophobic_annulus_3_aminos, hydrophobic_annulus_4_aminos]

## Set the parameters for each amino
for c in collection:
  for a in c: a.locate(test_file)

## Open files, where they exist, for structures
for index in range(equilibrated_point,end_frame):
  if(not isfile("FILE1/Stage1/s{}".format(index))): continue
  with open("FILE1/Stage1/s{}".format(index),"r") as f:
    flines = f.readlines()
    for c in collection:
      for a in c: a.add_data(flines)

## Calculate everything
f1 = open("Test_Strings","w")
f2 = open("Test_Atoms","w")
f3 = open("Test_Mus","w")
f4 = open("Test_Us","w")
f5 = open("Test_Traj","w")
for c in collection:
  for a in c: a.calculate_essentials(f1,f2,f3,f4,f5)

exit()

## Generates coefficients between aminos
generate_coefficient_dictionary(collection[0][0],collection[0][1])




with open("output_mus.csv","w") as f:
  for c in collection:
    for a in c:
      for i,j in zip(a.atomnames, a.mus):
        f.write("{},{},{},{},{}\n".format(a.string,i,j[0],j[1],j[2]))

with open("output_vs.csv","w") as f:
  for c in collection:
    for a in c:
      for i,j in zip(a.atomnames, a.v1s):
        f.write("{},{},{},{},{}\n".format(a.string,i,j[0],j[1],j[2]))


with open("output_corr.csv","w") as f:
  ## Calculate amino by amino correlation matrices
  for c1 in collection:
    for c2 in collection:
      #figure_index = 1
      #fig = plt.figure()
      for a1 in c1:
        for a2 in c2:
          ## Correlation matrices for each direction
          corr_mat11 = array([[ pearsonr(i,j)[0] for j in a2.projections1] for i in a1.projections1]) 
          corr_mat12 = array([[ pearsonr(i,j)[0] for j in a2.projections2] for i in a1.projections1]) 
          corr_mat13 = array([[ pearsonr(i,j)[0] for j in a2.projections3] for i in a1.projections1]) 
          corr_mat21 = array([[ pearsonr(i,j)[0] for j in a2.projections1] for i in a1.projections2]) 
          corr_mat22 = array([[ pearsonr(i,j)[0] for j in a2.projections2] for i in a1.projections2]) 
          corr_mat23 = array([[ pearsonr(i,j)[0] for j in a2.projections3] for i in a1.projections2]) 
          corr_mat31 = array([[ pearsonr(i,j)[0] for j in a2.projections1] for i in a1.projections3]) 
          corr_mat32 = array([[ pearsonr(i,j)[0] for j in a2.projections2] for i in a1.projections3]) 
          corr_mat33 = array([[ pearsonr(i,j)[0] for j in a2.projections3] for i in a1.projections3]) 
          for i1 in range(a1.numatoms):
            for i2 in range(a2.numatoms):
              f.write("{},{},{},{},{},11\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat11[i1,i2]))
              f.write("{},{},{},{},{},12\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat12[i1,i2]))
              f.write("{},{},{},{},{},13\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat13[i1,i2]))
              f.write("{},{},{},{},{},21\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat21[i1,i2]))
              f.write("{},{},{},{},{},22\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat22[i1,i2]))
              f.write("{},{},{},{},{},23\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat23[i1,i2]))
              f.write("{},{},{},{},{},31\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat31[i1,i2]))
              f.write("{},{},{},{},{},32\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat32[i1,i2]))
              f.write("{},{},{},{},{},33\n".format(a1.string,a1.atomnames[i1],a2.string,a2.atomnames[i2],corr_mat33[i1,i2]))
          #ax = fig.add_subplot(len(c1),len(c2),figure_index)
          #cax = ax.matshow(corr_mat)
          #fig.colorbar(cax)
          #if( figure_index < len(c2)): ax.set_title(a2.string)
          #ax.set_xticklabels([''] + a2.atomnames)
          #ax.set_yticklabels([''] + a1.atomnames)
          #ax.set_xlabel(a2.string, va='top')
          #if(figure_index % len(c2) == 1): ax.set_ylabel(a1.string)
          #figure_index += 1
      #plt.show()

## Calculate a list of coefficients K_{ab} for all atoms a and b
## Each coefficient is a matrix K

## This is natural to store in a double dictionary


global_dict = {}
for c1 in collection:
  for c2 in collection:
    for a1 in c1:
      for a2 in c2:
        ## A function to generate a dictionary of atom pairs
        global_dict.update(generate_coefficient_dictionary(a1,a2))
        



### Between each pair of atoms, calculate the coefficients



### Output a linear model

#We would like to output a linear model of the entire system under the impulses from a starting structure
#We solve this in the following way:

#we input the starting structure (mu) in x,y,z coordinates.
#we write the final position of each atom as mu_ij + dx_ij

#we calculate the dx_ij as functions of the impulse f_kl at all of the atoms!
#The coefficients are as follows

#Impulses are either in terms of (x,y,z) coordinates i.e. vectors at an atom or in terms of the principle motion of direction
#for the given atom.

#We convert the impulse at atom i, f_i = (fx,fy,fz) to a local impulse at the atom using the eigenvectors matrix V_i for the atom i as (V_i).f = (gx, gy, gz) = g_i

#We have a correlation matrix between principle directions in atom i and atom j, these alements are C_ij_xx, C_ij_xy, ... C_ij_z.
#We use this conversion to convert the impulse in local coordinates to a respose in local coordinates, the respose r is given by 
#rx = gx*C_ij_xx + gy*C_ij_yx + gz*C_ij_zx
#...
#rz = ...
#We need to sum this over all atoms,
#We need to convert back into x,y,z coordinates using the eigenvectors of the output atom!
#this gives us the displacement as a function of the impulse.










