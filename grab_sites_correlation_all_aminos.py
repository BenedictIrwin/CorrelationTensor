from numpy import float32, vectorize, array, transpose, corrcoef, mean, dot, matmul, subtract, save, load
from numpy.linalg import eigh, inv
from os.path import isfile
from scipy.stats.stats import pearsonr
#import scipy

import sys

## To convert a standard pdb string to a vector of x y z coords
def pdb_to_vec(s): return [float(s[30:38]),float(s[38:46]),float(s[46:54])]
def pdb_to_v(s): return [s[30:38].strip(),s[38:46].strip(),s[46:54].strip()]
def pdb_to_k(s): return s[17:26]
def to_vec(s): return [pdb_to_vec(i) for i in s]
def average_dic(x): return { i : mean(array(x[i]).astype(float),axis = 0) for i in x.keys()}


## Check the runtime input
argc = len(sys.argv)
if(argc!=3):
  print("Usage: python grab_sites...py [Equilibrated Frame e.g. 0] [End Frame e.g. 960]")
  exit()

## Make a non-class based version of the code

## Define a complete file that represents the structure of the protein for all frames
test_file = "Structures/s0"
equilibrated_point = int(sys.argv[1])
end_frame = int(sys.argv[2])

all_protein_aminos = []
num_atoms_dict = {}
previous_amino_key = "XYZABC"
count = 0
with open(test_file) as f:
  for line in f:
    if(line[0:3] == "TER"): continue
    if(line[0:3] == "END"): continue
    amino_key = line[17:26]
    if(amino_key==""): continue
    if(amino_key != previous_amino_key):
      count += 1
      previous_amino_key = amino_key
      all_protein_aminos.append(amino_key)
      num_atoms_dict[amino_key] = count
      count = 0
    else:
      count += 1

## Make sure files are valid
if isfile("valid_f.npy"):
  ## Allocate the required memory
  valid_files = load("valid_f.npy")
else:
  valid_files = []
  print("Checking structures, [this will only need to be run once!]")
  for index in range(equilibrated_point,end_frame):
    if(not isfile("Structures/s{}".format(index))): continue
    valid_files.append(index)
  save("valid_f",valid_files)
  valid_files = load("valid_f.npy")
  print(valid_files)


data = []

## This list of terminal sites is protein specific!
terminal_sites = ["","HIS A 306","GLU A 360","ASP A 665","THR A 720","LYS A1027","MET A1093","HIS A1398","THR A1452","LYS A1758","GLN A1823"]



print("Number of valid files: ",len(valid_files))

if isfile("Trajectory.npy"):
  bbb = load("Trajectory.npy")
  bbb = transpose(bbb, (1,2,0))

else:
  print("Generating Trajectory: This may take quite some time for large simulations... but will only need to be done once.") 
  for index in valid_files:
    dic = { i : [] for i in all_protein_aminos}
    with open("Structures/s{}".format(index)) as f:
      flines = f.readlines()
      xs = [ pdb_to_v(i) for i in flines]
      keys = [ pdb_to_k(i) for i in flines]
    for i,j in zip(keys,xs):
      if(i in terminal_sites): continue
      dic[i].append(j)
  
    data.append(list(average_dic(dic).values()))
    bbb = array(data)
    save("Trajectory", bbb)
    bbb = transpose(bbb, (1,2,0))


#save("BIG_Labels",all_protein_aminos)
print("Trajectory Shape: ",bbb.shape)
#aaa = load("BIG_labels.npy")
## Shape = 75000, 1823, 3 

num_aminos = bbb.shape[0]
print("Num Aminos: ", num_aminos)

if(True):
  print("Generating Means")
  mus = []
  for i in bbb:
    mus.append( mean(i, axis = 1))
  mus=array(mus)
  print(mus)
  print(mus.shape)
  save("Means",array(mus))

if(True):
  print("Generating Spatial Correlation Matrices")
  Us=[]
  Uinvs = []
  for i in bbb:
    try:
      _, v = eigh(corrcoef(i))
    except:
      print(corrcoef(i))
    #  _, v = scipy.linalg.eigh(corrcoef(i))
    ## The eigenvalue matrix used to project
    Us.append(v[:,::-1])
    Uinvs.append(transpose(v[:,::-1]))
  Us = array(Us)
  Uinvs = array(Uinvs)
  print(Us)
  print(Us.shape) 
  save("SpatialCorrelation",array(Us))

  print("Generating Time Projections")
  Trajs = []
  for i in range(num_aminos):
    d = transpose(bbb[i])
    mu = mus[i]
    qqq = matmul(Uinvs[i],transpose(subtract(d,mu))) 
    Trajs.append(qqq)

  Trajs = array(Trajs)
  save("Projections",array(Trajs))
  

exit()

save("GIANTModel_Strings",array(Strings))
save("GIANTModel_Atoms",array(Atoms))

#    ## Calculate projections
#    ## Each atom has three projections
#    ## Two atoms have 3x3 combinations of projections
#    for i in range(self.num_atoms):
#      d = transpose(self.store[i])
#      mu = self.mus[i]
#      qqq = matmul(self.Uinvs[i],transpose(subtract(d,mu)))
#      name = self.atomnames[i]
#
#      S.append(self.string)
#      A.append(name)
#      M.append(mu)
#      U.append(self.Us[i])
#      T.append(qqq)


exit()

## Define a list as a global collection
#collection = [binding_site_aminos, top_annulus_aminos, middle_annulus_aminos, bottom_annulus_aminos, OH_annulus_1_aminos, OH_annulus_2_aminos, OH_annulus_3_aminos, OH_annulus_4_aminos, hydrophobic_annulus_1_aminos, hydrophobic_annulus_2_aminos, hydrophobic_annulus_3_aminos, hydrophobic_annulus_4_aminos]

collection = [[Amino(s) for s in all_protein_aminos]]

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

Strings = []
Atoms = []
Mus = []
Us = []
Trajs = []
for c in collection:
  for a in c: a.calculate_essentials(Strings,Atoms,Mus,Us,Trajs)




NN = 12852 

#pdb_vec = vectorize(pdb_to_vec)

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
    self.store.append( to_vec(flines[self.start:self.stop]) )

  ## Calculate the principle direction of motion
  ## Calculate the size of the eigenvalue
  ## Calculate the projections
  def calculate_essentials(self,S,A,M,U,T):

    self.store = transpose( array(self.store), (1,2,0) )
    self.num_atoms = self.store.shape[0]
    for i in self.store:
      self.mus.append(mean(i, axis = 1))
      ## Hermitian matrix
      try:
        _, v = eigh(corrcoef(i))
      except:
        print(corrcoef(i))
      #  _, v = scipy.linalg.eigh(corrcoef(i))
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

      S.append(self.string)
      A.append(name)
      M.append(mu)
      U.append(self.Us[i])
      T.append(qqq)

 


