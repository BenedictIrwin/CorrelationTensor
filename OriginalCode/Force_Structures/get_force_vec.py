import numpy as np

### apo --> holo is adding GABA
### holo --> apo is removing GABA

apo_file = "apo_GABAAR_fitted.pdb"
holo_file = "holo_GABAAR_fitted.pdb"
s_file = "s25000"

## Get the common set of aminos between the two files

apo_file_dic = {}
with open(apo_file,"r") as f:
  for line in f:
    line = line.strip()
    amino = line[17:26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    arr = np.array([1,x,y,z])    
    if(amino in apo_file_dic.keys()): apo_file_dic[amino]+= arr
    else: apo_file_dic[amino] = arr
## Average over numebr of atoms
apo_file_dic = {i : apo_file_dic[i][1:]/apo_file_dic[i][0] for i in apo_file_dic.keys()}

holo_file_dic = {}
with open(holo_file,"r") as f:
  for line in f:
    line = line.strip()
    amino = line[17:26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    arr = np.array([1,x,y,z])    
    if(amino in holo_file_dic.keys()): holo_file_dic[amino]+= arr
    else: holo_file_dic[amino] = arr
## Average over numebr of atoms
holo_file_dic = {i : holo_file_dic[i][1:]/holo_file_dic[i][0] for i in holo_file_dic.keys()}


apo_keys_set = set(apo_file_dic.keys())
holo_keys_set = set(holo_file_dic.keys())

## Common keys 

common_keys = apo_keys_set.intersection(holo_keys_set)

## Generate forward (binding differences)
forward_force_dic = {}
backward_force_dic = {}
for i in common_keys:
   forward_force_dic[i] = holo_file_dic[i] - apo_file_dic[i] 
   backward_force_dic[i] = apo_file_dic[i] - holo_file_dic[i] 


average_forward_shift = np.mean([ i for i in forward_force_dic.values()], axis = 0)
average_backward_shift = np.mean([ i for i in backward_force_dic.values()], axis = 0)

## Shifts
print(average_forward_shift)
print(average_backward_shift)

## Optional coordinate sets
corrected_forward_force_dic = {i : forward_force_dic[i] - average_forward_shift for i in forward_force_dic.keys()}
corrected_backward_force_dic = {i : backward_force_dic[i] - average_backward_shift for i in backward_force_dic.keys()}


## Get the keys to map ontop s25000

convert_s_to_a = {"PHE I1504":"PHE A  65","ARG I1506":"ARG A  67","LEU I1557":"LEU A 118", "THR I1569":"THR A 130","PHE E 773":"PHE D  65","ARG E 775":"ARG D  67","LEU E 826":"LEU D 118", "THR E 838":"THR D 130", "TYR A  89":"TYR B  97","GLU A 147":"GLU B 155", "SER A 148":"SER B 156", "TYR A 149":"TYR B 157", "PHE A 192":"PHE B 200", "SER A 193":"ALA B 201", "THR A 194":"THR B 202", "TYR A 197":"TYR B 205", "TYR G1181":"TYR E  97","GLU G1239":"GLU E 155", "SER G1240":"SER E 156", "TYR G1241":"TYR E 157", "PHE G1284":"PHE E 200", "SER G1285":"ALA E 201", "THR G1286":"THR E 202", "TYR G1289":"TYR E 205"}

## The bigger 6AA dictionary
convert_s_to_a = {"PHE I1504":"PHE A  65","PHE I1505":"PHE A  66","ARG I1506":"ARG A  67","LEU I1557":"LEU A 118","THR I1569":"THR A 130","MET I1570":"MET A 131","PHE E 773":"PHE D  65","PHE E 774":"PHE D  66","ARG E 775":"ARG D  67","LEU E 826":"LEU D 118","THR E 838":"THR D 130","MET E 839":"MET D 131","TYR A  89":"TYR B  97","GLU A 147":"GLU B 155","SER A 148":"SER B 156","TYR A 149":"TYR B 157","GLY A 150":"GLY B 158","PHE A 192":"PHE B 200","SER A 193":"ALA B 201","THR A 194":"THR B 202","GLY A 195":"GLY B 203","TYR A 197":"TYR B 205","TYR G1181":"TYR E  97","GLU G1239":"GLU E 155","SER G1240":"SER E 156","TYR G1241":"TYR E 157","GLY G1242":"GLY E 158","PHE G1284":"PHE E 200","SER G1285":"ALA E 201","THR G1286":"THR E 202","GLY G1287":"GLY E 203","TYR G1289":"TYR E 205"}

## Inverse dictionary
convert_a_to_s = { convert_s_to_a[i] : i for i in convert_s_to_a.keys()}

## Sites 1 and 2 for GABA
## Define GABA site 1 as A,B+I,J which is A and B in Aricescu
## Define GABA site 2 as E,F+G,H which is D and E in Aricescu
binding_site_aminos = ["PHE A  65", "ARG A  67", "PHE D  65", "ARG D  67", "TYR B 157", "TYR B 205", "TYR E 157", "TYR E 205"]
binding_site_2_aminos = ["PHE D  65", "ARG D  67","TYR E 157", "TYR E 205"]
binding_site_1_aminos = ["PHE A  65", "ARG A  67","TYR B 157", "TYR B 205"]
#### 5AA acids
AA5_acids = ["PHE A  65", "ARG A  67","LEU A 118","THR A 130","TYR B  97","GLU B 155","SER B 156","TYR B 157","PHE B 200","ALA B 201","THR B 202","TYR B 205","PHE D  65","ARG D  67","LEU D 118","THR D 130","TYR E  97","GLU E 155","SER E 156","TYR E 157","PHE E 200","ALA E 201","THR E 202","TYR E 205"]

AA5_site_1_acids = ["PHE A  65", "ARG A  67","LEU A 118","THR A 130","TYR B  97","GLU B 155","SER B 156","TYR B 157","PHE B 200","ALA B 201","THR B 202","TYR B 205"]
AA5_site_2_acids = ["PHE D  65","ARG D  67","LEU D 118","THR D 130","TYR E  97","GLU E 155","SER E 156","TYR E 157","PHE E 200","ALA E 201","THR E 202","TYR E 205"]

#### 4AA acids
AA4_acids = ["PHE A  65", "ARG A  67","THR A 130","TYR B  97","GLU B 155","SER B 156","TYR B 157","THR B 202","TYR B 205","PHE D  65","ARG D  67","THR D 130","TYR E  97","GLU E 155","SER E 156","TYR E 157","THR E 202","TYR E 205"]
AA4_site_1_acids = ["PHE A  65", "ARG A  67","THR A 130","TYR B  97","GLU B 155","SER B 156","TYR B 157","THR B 202","TYR B 205"]
AA4_site_2_acids = ["PHE D  65","ARG D  67","THR D 130","TYR E  97","GLU E 155","SER E 156","TYR E 157","THR E 202","TYR E 205"]


AA6_acids = ["PHE A  65","PHE A  66", "ARG A  67","LEU A 118","THR A 130","MET A 131","PHE D  65","PHE D  66","ARG D  67","LEU D 118","THR D 130","MET D 131","TYR B  97","GLU B 155","SER B 156","TYR B 157","GLY B 158","PHE B 200","ALA B 201","THR B 202","GLY B 203","TYR B 205","TYR E  97","GLU E 155","SER E 156","TYR E 157","GLY E 158","PHE E 200","ALA E 201","THR E 202","GLY E 203","TYR E 205"]
AA6_site_1_acids = ["PHE A  65","PHE A  66", "ARG A  67","LEU A 118","THR A 130","MET A 131","TYR B  97","GLU B 155","SER B 156","TYR B 157","GLY B 158","PHE B 200","ALA B 201","THR B 202","GLY B 203","TYR B 205"]
AA6_site_2_acids = ["PHE D  65","PHE D  66","ARG D  67","LEU D 118","THR D 130","MET D 131","TYR E  97","GLU E 155","SER E 156","TYR E 157","GLY E 158","PHE E 200","ALA E 201","THR E 202","GLY E 203","TYR E 205"]


## 9 groups, minimal: 1,2,12, 4AA: 1,2,12, 5AA: 1,2,12
## 4 directions: forward, backward, forward_corrected, backward_corrected
## This is 36 calculations

binding_site_on_s = np.array([convert_a_to_s[i] for i in binding_site_aminos])
binding_site_1_on_s = [convert_a_to_s[i] for i in binding_site_1_aminos]
binding_site_2_on_s = [convert_a_to_s[i] for i in binding_site_2_aminos]

AA4_acids_on_s = np.array([convert_a_to_s[i] for i in AA4_acids])
AA4_site_1_acids_on_s = np.array([convert_a_to_s[i] for i in AA4_site_1_acids])
AA4_site_2_acids_on_s = np.array([convert_a_to_s[i] for i in AA4_site_2_acids])

AA5_acids_on_s = [convert_a_to_s[i] for i in AA5_acids]
AA5_site_1_acids_on_s = [convert_a_to_s[i] for i in AA5_site_1_acids]
AA5_site_2_acids_on_s = [convert_a_to_s[i] for i in AA5_site_2_acids]

AA6_acids_on_s = [convert_a_to_s[i] for i in AA6_acids]
AA6_site_1_acids_on_s = [convert_a_to_s[i] for i in AA6_site_1_acids]
AA6_site_2_acids_on_s = [convert_a_to_s[i] for i in AA6_site_2_acids]


s_strings = np.load("../BIG_Labels.npy")


## Minimal
I_S12_min = [ np.argwhere(s_strings == i)[0][0] for i in binding_site_on_s]
F_S12_min_f = [ forward_force_dic[i] for i in binding_site_aminos ]
#F_S12_min_fc = [ corrected_forward_force_dic[i] for i in binding_site_aminos ]
F_S12_min_b = [ backward_force_dic[i] for i in binding_site_aminos ]
#F_S12_min_bc = [ corrected_backward_force_dic[i] for i in binding_site_aminos ]

I_S1_min = [ np.argwhere(s_strings == i)[0][0] for i in binding_site_1_on_s]
F_S1_min_f = [ forward_force_dic[i] for i in binding_site_1_aminos ]
#F_S1_min_fc = [ corrected_forward_force_dic[i] for i in binding_site_1_aminos ]
F_S1_min_b = [ backward_force_dic[i] for i in binding_site_1_aminos ]
#F_S1_min_bc = [ corrected_backward_force_dic[i] for i in binding_site_1_aminos ]

I_S2_min = [ np.argwhere(s_strings == i)[0][0] for i in binding_site_2_on_s]
F_S2_min_f = [ forward_force_dic[i] for i in binding_site_2_aminos ]
#F_S2_min_fc = [ corrected_forward_force_dic[i] for i in binding_site_2_aminos ]
F_S2_min_b = [ backward_force_dic[i] for i in binding_site_2_aminos ]
#F_S2_min_bc = [ corrected_backward_force_dic[i] for i in binding_site_2_aminos ]

## 4AA
I_S12_4AA = [ np.argwhere(s_strings == i)[0][0] for i in AA4_acids_on_s]
F_S12_4AA_f = [ forward_force_dic[i] for i in AA4_acids ]
#F_S12_4AA_fc = [ corrected_forward_force_dic[i] for i in AA4_acids ]
F_S12_4AA_b = [ backward_force_dic[i] for i in AA4_acids ]
#F_S12_4AA_bc = [ corrected_backward_force_dic[i] for i in AA4_acids ]

I_S1_4AA = [ np.argwhere(s_strings == i)[0][0] for i in AA4_site_1_acids_on_s]
F_S1_4AA_f = [ forward_force_dic[i] for i in AA4_site_1_acids ]
#F_S1_4AA_fc = [ corrected_forward_force_dic[i] for i in AA4_site_1_acids ]
F_S1_4AA_b = [ backward_force_dic[i] for i in AA4_site_1_acids ]
#F_S1_4AA_bc = [ corrected_backward_force_dic[i] for i in AA4_site_1_acids ]

I_S2_4AA = [ np.argwhere(s_strings == i)[0][0] for i in AA4_site_2_acids_on_s]
F_S2_4AA_f = [ forward_force_dic[i] for i in AA4_site_2_acids ]
#F_S2_4AA_fc = [ corrected_forward_force_dic[i] for i in AA4_site_2_acids ]
F_S2_4AA_b = [ backward_force_dic[i] for i in AA4_site_2_acids ]
#F_S2_4AA_bc = [ corrected_backward_force_dic[i] for i in AA4_site_2_acids ]

## 5AA
I_S12_5AA = [ np.argwhere(s_strings == i)[0][0] for i in AA5_acids_on_s]
F_S12_5AA_f = [ forward_force_dic[i] for i in AA5_acids ]
#F_S12_5AA_fc = [ corrected_forward_force_dic[i] for i in AA5_acids ]
F_S12_5AA_b = [ backward_force_dic[i] for i in AA5_acids ]
#F_S12_5AA_bc = [ corrected_backward_force_dic[i] for i in AA5_acids ]

I_S1_5AA = [ np.argwhere(s_strings == i)[0][0] for i in AA5_site_1_acids_on_s]
F_S1_5AA_f = [ forward_force_dic[i] for i in AA5_site_1_acids ]
#F_S1_5AA_fc = [ corrected_forward_force_dic[i] for i in AA5_site_1_acids ]
F_S1_5AA_b = [ backward_force_dic[i] for i in AA5_site_1_acids ]
#F_S1_5AA_bc = [ corrected_backward_force_dic[i] for i in AA5_site_1_acids ]

I_S2_5AA = [ np.argwhere(s_strings == i)[0][0] for i in AA5_site_2_acids_on_s]
F_S2_5AA_f = [ forward_force_dic[i] for i in AA5_site_2_acids ]
#F_S2_5AA_fc = [ corrected_forward_force_dic[i] for i in AA5_site_2_acids ]
F_S2_5AA_b = [ backward_force_dic[i] for i in AA5_site_2_acids ]
#F_S2_5AA_bc = [ corrected_backward_force_dic[i] for i in AA5_site_2_acids ]

## 6AA
I_S12_6AA = [ np.argwhere(s_strings == i) for i in AA6_acids_on_s]

print(I_S12_6AA)
print(AA6_acids_on_s)

I_S12_6AA = [ np.argwhere(s_strings == i)[0][0] for i in AA6_acids_on_s]
F_S12_6AA_f = [ forward_force_dic[i] for i in AA6_acids ]
F_S12_6AA_b = [ backward_force_dic[i] for i in AA6_acids ]

I_S1_6AA = [ np.argwhere(s_strings == i)[0][0] for i in AA6_site_1_acids_on_s]
F_S1_6AA_f = [ forward_force_dic[i] for i in AA6_site_1_acids ]
F_S1_6AA_b = [ backward_force_dic[i] for i in AA6_site_1_acids ]

I_S2_6AA = [ np.argwhere(s_strings == i)[0][0] for i in AA6_site_2_acids_on_s]
F_S2_6AA_f = [ forward_force_dic[i] for i in AA6_site_2_acids ]
F_S2_6AA_b = [ backward_force_dic[i] for i in AA6_site_2_acids ]

## Make sure to print out the shifts!
## Make sure to do a p-test of binding site vs the entire distribution!

if(False):
  print("Shifts")
  for i,j,k in zip(F_S12_5AA_f,I_S12_5AA_f,F_S12_5AA_fc):
    print("{0:s},{1:.3f},{2:.3f}".format(s_strings[j],np.sqrt(np.dot(i,i)),np.sqrt(np.dot(k,k))))
  
  from matplotlib import pyplot as plt
  
  ### Make a histogram of force dic values:
  all_values = [np.sqrt(np.dot(i,i)) for i in forward_force_dic.values()]
  plt.hist(all_values,bins =20)
  plt.xlabel("Shift in Angstroms")
  plt.ylabel("Frequency")
  plt.title("Forward Shift")
  plt.show()
  
  all_values = [np.sqrt(np.dot(i,i)) for i in corrected_forward_force_dic.values()]
  plt.hist(all_values,bins =20)
  plt.xlabel("Shift in Angstroms")
  plt.ylabel("Frequency")
  plt.title("Forward Shift Corrected")
  plt.show()

n_amino = len(s_strings)

def save_force(F,I,tag):
  forces = np.zeros(shape = (n_amino,3))
  for force,index in zip(F,I):
    forces[index]=force
  np.save("GIANTModel_Forces_{}".format(tag),forces)

## Save all of the forces in files

### Min
save_force(F_S12_min_f,I_S12_min,"min_f12")
#save_force(F_S12_min_fc,I_S12_min,"min_fc12")
save_force(F_S12_min_b,I_S12_min,"min_b12")
#save_force(F_S12_min_bc,I_S12_min,"min_bc12")

save_force(F_S1_min_f,I_S1_min,"min_f1")
#save_force(F_S1_min_fc,I_S1_min,"min_fc1")
save_force(F_S1_min_b,I_S1_min,"min_b1")
#save_force(F_S1_min_bc,I_S1_min,"min_bc1")

save_force(F_S2_min_f,I_S2_min,"min_f2")
#save_force(F_S2_min_fc,I_S2_min,"min_fc2")
save_force(F_S2_min_b,I_S2_min,"min_b2")
#save_force(F_S2_min_bc,I_S2_min,"min_bc2")


save_force(F_S12_4AA_f,I_S12_4AA,"4AA_f12")
#save_force(F_S12_4AA_fc,I_S12_4AA,"4AA_fc12")
save_force(F_S12_4AA_b,I_S12_4AA,"4AA_b12")
#save_force(F_S12_4AA_bc,I_S12_4AA,"4AA_bc12")

save_force(F_S1_4AA_f,I_S1_4AA,"4AA_f1")
#save_force(F_S1_4AA_fc,I_S1_4AA,"4AA_fc1")
save_force(F_S1_4AA_b,I_S1_4AA,"4AA_b1")
#save_force(F_S1_4AA_bc,I_S1_4AA,"4AA_bc1")

save_force(F_S2_4AA_f,I_S2_4AA,"4AA_f2")
#save_force(F_S2_4AA_fc,I_S2_4AA,"4AA_fc2")
save_force(F_S2_4AA_b,I_S2_4AA,"4AA_b2")
#save_force(F_S2_4AA_bc,I_S2_4AA,"4AA_bc2")


save_force(F_S12_5AA_f,I_S12_5AA,"5AA_f12")
#save_force(F_S12_5AA_fc,I_S12_5AA,"5AA_fc12")
save_force(F_S12_5AA_b,I_S12_5AA,"5AA_b12")
#save_force(F_S12_5AA_bc,I_S12_5AA,"5AA_bc12")

save_force(F_S1_5AA_f,I_S1_5AA,"5AA_f1")
#save_force(F_S1_5AA_fc,I_S1_5AA,"5AA_fc1")
save_force(F_S1_5AA_b,I_S1_5AA,"5AA_b1")
#save_force(F_S1_5AA_bc,I_S1_5AA,"5AA_bc1")

save_force(F_S2_5AA_f,I_S2_5AA,"5AA_f2")
#save_force(F_S2_5AA_fc,I_S2_5AA,"5AA_fc2")
save_force(F_S2_5AA_b,I_S2_5AA,"5AA_b2")
#save_force(F_S2_5AA_bc,I_S2_5AA,"5AA_bc2")


## Savethe 6AA forces
save_force(F_S12_6AA_f,I_S12_6AA,"6AA_f12")
save_force(F_S12_6AA_b,I_S12_6AA,"6AA_b12")

save_force(F_S1_6AA_f,I_S1_6AA,"6AA_f1")
save_force(F_S1_6AA_b,I_S1_6AA,"6AA_b1")

save_force(F_S2_6AA_f,I_S2_6AA,"6AA_f2")
save_force(F_S2_6AA_b,I_S2_6AA,"6AA_b2")



