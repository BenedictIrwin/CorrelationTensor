#!/bin/python

import sys
argc=len(sys.argv)

filename=sys.argv[1]
target="END"

outfix="s"
i=0

g=open(outfix+str(i),"w")
with open(filename) as f:
  for fullline in f:
    line=fullline.split()
    if(line[0]=="ATOM" or line[0]=="TER"): print(fullline,file=g,end="")
    if(line[0]==target):
      print(fullline,file=g,end="")
      g.close()
      i+=1
      g=open(outfix+str(i),"w")
g.close()
