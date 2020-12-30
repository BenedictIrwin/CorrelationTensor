import sys
import numpy as np
import os
## A script to assemble a set of .pdb frames
## The frames represent the correlated motion under a custom input force

## Writes a PDB frame with only the MU positions
#def mus_to_pdb(fn,m,a,s,vs,sc):
#  with open(fn,"w") as f:
#    kk = 0
#    for mus,aa,ss,vv in zip(m,a,s,vs):
#      x = mus + sc*vv 
#      f.write("ATOM  {0: >5d}  {2:}  {1:}    {3:8.3f}{4:8.3f}{5:8.3f}\n".format(kk,ss,aa,x[0],x[1],x[2]))


argc = len(sys.argv)
if( argc != 2):
  print("Usage : [A name for these files]")
  exit()

file_name = sys.argv[1]
tag = file_name.replace("GIANTModel_Deltas_","").replace(".npy","")
print(tag)
out_name = tag

input_pdb = "s25000"
Deltas = np.load(file_name)
labels = np.load("BIG_Labels.npy")
amino_dic = {labels[i] : i for i in range(len(labels))}
shifts = [np.round(0.1 * i,2) for i in range(10+1)]

for sh in shifts:
  os.system("rm {}_{}.pdb".format(out_name,sh))


olines = []
ss = []
origs = []
with open(input_pdb,"r") as f:
  for line in f:
    line = line.strip()
    if(line[0:3]=="TER" or line[17:26] == ""): continue
    olines.append(line[:30])
    ss.append(Deltas[amino_dic[line[17:26]]])
    origs.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))


#if line[0:3] == "TER" : continue
#oline = line[:30]
#amino = line[17:26]
#if amino == "" : continue
#shift = Deltas[amino_dic[amino]]
#orig = np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])


for sh in shifts:
  print(sh)
  with open("{}_{}.pdb".format(out_name,sh),"a") as g:
    for o,s,og in zip(olines,ss,origs):
      new = og + sh*s
      g.write(o+"{0:8.3f}{1:8.3f}{2:8.3f}\n".format(new[0],new[1],new[2]))

os.system("rm {}_cat.pdb".format(out_name))
for sh in shifts:
  os.system("cat {}_{}.pdb >> {}_cat.pdb".format(out_name,sh,out_name))
  os.system("echo \'ENDMDL\' >> {}_cat.pdb".format(out_name))


