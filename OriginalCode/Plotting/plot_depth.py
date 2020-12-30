import numpy as np
import sys
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

tag = sys.argv[1]

indices = ["0.0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"]
files = [ tag+"_"+j+".pdbdepth" for j in indices ]


data= []

for file_name in files:
  with open(file_name,"r") as f: lines = f.readlines()
  lines = np.array([ i.strip().split(" ") for i in lines  ])
  lines = np.transpose(lines).astype(float)
  lines = np.array(sorted(lines, key = lambda x: x[0]))
  data.append(lines[1])

data = np.transpose(np.array(data))


opens = []
closes = []
highs = []
lows = []
for i in data:
  opens.append(i[0])
  closes.append(i[-1])
  highs.append(max(i))
  lows.append(min(i))

opens = np.array(opens)

from matplotlib.finance import candlestick2_ochl



tag_to_dist = {"T" : 19.95045714285714,"M" : 13.535514285714285,"B" : -22.043857142857142,"OH1" : 15.330485714285714,"OH2" : 1.617575,"OH3" : -2.5823,"OH4" : -8.916849999999998,"HP1" : 17.293228571428568,"HP2" : 3.939978571428571,"HP3" : -4.974257142857143,"HP4" : -12.003725,"CC" : -1.9202571428571429}



## Code for the absolute
if(True):
  bc = 'white'
  fc = 'black'

  fig = plt.figure(facecolor = bc)
  fig.set_size_inches(18.5, 10.5)
  ax1 = plt.subplot2grid((1,1), (0,0))
  ax1.set_facecolor(bc) 

  ax1.spines['bottom'].set_color(fc)
  ax1.spines['top'].set_color(fc) 
  ax1.spines['right'].set_color(fc)
  ax1.spines['left'].set_color(fc)

  ax1.tick_params(axis='x', colors=fc)
  ax1.tick_params(axis='y', colors=fc)

 
  candlestick2_ochl(ax1,opens,closes,highs,lows,width= 1, colorup='#77d879', colordown='#db3f3f')
  plt.xlabel('Z-coordinate [Angstrom]', fontsize = 24, color = fc)
  plt.xticks([i for i in range(len(lines[0]))][::5],lines[0][::5], fontsize = 18)
  plt.yticks(fontsize=18)
  plt.ylabel('Smoothed Channel Radius [Angstrom]',fontsize = 24, color = fc)
  
  titletag = tag.replace("min_","Minimal").replace("6AA_","Residues within 6 $\AA$").replace("4AA_","Residues within 4 $\AA$").replace("5AA_","Residues within 5 $\AA$").replace("fc"," Forward (Correction)").replace("bc", " Backward (Correction)").replace("f"," Forward").replace("b"," Backward").replace("12"," Sites AB and DE").replace("1"," Site AB").replace("2"," Site DE")
  plt.title(titletag, fontsize = 24, color = fc)
  plt.subplots_adjust(left=0.09, bottom=0.20, right=0.94, top=0.90, wspace=0.2, hspace=0)
 
   
  for i,j in zip(tag_to_dist.keys(),tag_to_dist.values()):
    plt.plot([j+40, j+40], [1, 9], 'k--', lw=1, color =fc)
    shif = 0.0
    if(i=="CC"): shif = 0.4
    if(i=="OH3"): shif = -0.4
    plt.text(j+40-0.8+shif,9.05,i, fontsize = 18, rotation = 90, ha = "left", va = "bottom", color = fc)

  #print(tag_to_dist.values())
  #exit()
  new_x = np.linspace(0,100,300)
  open_spline = interp1d(range(0,101),opens,kind="cubic")
  close_spline = interp1d(range(0,101),closes,kind="cubic")
  plt.plot(new_x,open_spline(new_x),'-',color='black',linewidth=1,label="Original Pore")
  plt.plot(new_x,close_spline(new_x),'-',color='black',linewidth=2, label="New Pore")
  plt.legend(loc="lower right", fontsize=18)
  #plt.show()
  #exit() 
  fig.savefig('{}_LARGE.png'.format(tag), dpi=300)
  exit() 

## Code for the shifted
if(False):
  fig = plt.figure()
  fig.set_size_inches(18.5, 10.5)
  ax1 = plt.subplot2grid((1,1), (0,0))
  
  candlestick2_ochl(ax1,[ 0 for i in range(len(opens))],np.array(closes)-np.array(opens),np.array(highs)-np.array(opens),np.array(lows)-np.array(opens),width= 1, colorup='#77d879', colordown='#db3f3f')
  plt.xlabel('Z-coordinate [Angstrom]', fontsize = 18)
  plt.xticks([i for i in range(len(lines[0]))][::5],lines[0][::5], fontsize = 14)
  plt.ylabel('Smoothed Channel Radius [Angstrom]',fontsize = 18)
  
  titletag = tag.replace("min_","Minimal").replace("6AA_","Residues within 6 $\AA$").replace("4AA_","Residues within 4 $\AA$").replace("5AA_","Residues within 5 $\AA$").replace("fc"," Forward (Correction)").replace("bc", " Backward (Correction)").replace("f"," Forward").replace("b"," Backward").replace("12"," Sites AB and DE").replace("1"," Site AB").replace("2"," Site DE")
  plt.title(titletag, fontsize = 18)
  plt.legend()
  plt.subplots_adjust(left=0.09, bottom=0.20, right=0.94, top=0.90, wspace=0.2, hspace=0)
  #plt.show()
  for i,j in zip(tag_to_dist.keys(),tag_to_dist.values()):
    plt.plot([j+40, j+40], [-1, 0.9], 'k--', lw=1)
    shif = 0.0
    if(i=="CC"): shif = 0.25
    if(i=="OH3"): shif = -0.25
    plt.text(j+40-0.5+shif,0.94,i, fontsize = 14, rotation = 90, ha = "left", va = "bottom")
  fig.savefig('{}_relative.png'.format(tag), dpi=100)


