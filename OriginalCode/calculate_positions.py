from numpy import load, einsum, save
import sys

file_name = sys.argv[1]
tag = file_name.replace("GIANTModel_Forces_","").replace(".npy","")
print(tag)

## Read in the required varaibles
Forces = load(file_name)
Coeffs = load("GIANTModel_Coeffs.npy")

### Calculate
Deltas = einsum("ai,abij->bj",Forces,Coeffs)

## Save
save("GIANTModel_Deltas_{}".format(tag),Deltas)

