Crankite Data Processing
Benedict Irwin 2018

Stage 1:

The first stage is to take the crankite output which is one very long .pdb file and spit it into multiple individual .pdb files
These files will be used by the clustering algorithm.

Run the program "splitfile.py" on the large .pdb file, this will save many files called
s0
s1
s2
....
sN

where N is the number of snapshots in the crankite trajectory, somethign like 100,000
