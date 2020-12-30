import numpy as np
import sys
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

vert_shift = 0
bc = 'white'
fc = 'black'

#fig = plt.figure(facecolor = bc)
#fig.set_size_inches(18.5, 18.5)
#ax1 = plt.subplot2grid((1,1), (0,0))
#ax1.set_facecolor(bc) 

fig, axs = plt.subplots(6,1,sharex=True)

fig.set_size_inches(11.69,8.27)
fig.set_size_inches(8.27,11.69)
fig.set_size_inches(2*8.27,2*11.69)

#ax1.spines['bottom'].set_color(fc)
#ax1.spines['top'].set_color(fc) 
#ax1.spines['right'].set_color(fc)
#ax1.spines['left'].set_color(fc)

#ax1.tick_params(axis='x', colors=fc)
#ax1.tick_params(axis='y', colors=fc)
## Loop over tags
indd = 0
for tag, ax in zip(["min_f1","min_f2","min_f12","6AA_f1","6AA_f2","6AA_f12"],axs):
  
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
  
  opens = np.array(opens) - vert_shift
  closes = np.array(closes) - vert_shift
  highs = np.array(highs) - vert_shift
  lows = np.array(lows) - vert_shift
  
  from matplotlib.finance import candlestick2_ochl
  
  
  
  tag_to_dist = {"T" : 19.95045714285714,"M" : 13.535514285714285,"B" : -22.043857142857142,"OH1" : 15.330485714285714,"OH2" : 1.617575,"OH3" : -2.5823,"OH4" : -8.916849999999998,"HP1" : 17.293228571428568,"HP2" : 3.939978571428571,"HP3" : -4.974257142857143,"HP4" : -12.003725,"CC" : -1.9202571428571429}
  
  
  
  ## Code for the absolute
  if(True):
  
   
    candlestick2_ochl(ax,opens,closes,highs,lows,width= 1, colorup='#77d879', colordown='#db3f3f')
    if(tag=="6AA_f12"):
      ax.set_xlabel('Z-coordinate [Angstrom]', fontsize = 24, color = fc)
      #ax.set_xticks([i for i in range(len(lines[0]))][::10],lines[0][::10])
      ax.set_xticks(np.array([5*i for i in range(21)]))
      ax.set_xticklabels(np.array([5*i-40 for i in range(21)]))
      #ax.set_yticks(fontsize=22)
      ax.set_ylabel('.                                                        Smoothed Channel Radius [Angstrom]',fontsize = 24, color = fc)
    
    titletag = tag.replace("min_","Minimal").replace("6AA_","Residues within 6 $\AA$").replace("4AA_","Residues within 4 $\AA$").replace("5AA_","Residues within 5 $\AA$").replace("fc"," Forward (Correction)").replace("bc", " Backward (Correction)").replace("f"," Forward").replace("b"," Backward").replace("12"," Sites AB and DE").replace("1"," Site AB").replace("2"," Site DE")
    #plt.title(titletag, fontsize = 24, color = fc)
    #ax.subplots_adjust(left=0.09, bottom=0.20, right=0.94, top=0.90, wspace=0.2, hspace=0)
   
     
    #if(tag=="min_f1"):
    ax.set_yticks(np.array([2.5+2.5*i for i in range(4)]))
    ax.set_ylim((1.5,10.0))
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    #for i,j in zip(tag_to_dist.keys(),tag_to_dist.values()):
    #  ax.plot([j+40, j+40], [-29.5, 12], 'k--', lw=1, color =fc)
    #  shif = 0.0
    #  if(i=="CC"): shif = 0.4
    #  if(i=="OH3"): shif = -0.4
    #  if(tag=="min_f1"):
    #    ax.text(j+40-0.8+shif,10.5,i, fontsize = 18, rotation = 90, ha = "left", va = "bottom", color = fc)

    if(tag == "min_f1"): minilabel = "Min. AB" 
    if(tag == "min_f2"): minilabel = "Min. DE" 
    if(tag == "min_f12"): minilabel = "Min. AB+DE" 
    if(tag == "6AA_f1"): minilabel = "6$\AA$ AB" 
    if(tag == "6AA_f2"): minilabel = "6$\AA$ DE" 
    if(tag == "6AA_f12"): minilabel = "6$\AA$ AB+DE" 
    #ax.text(90,8.0,minilabel,fontsize=18, color = fc)
 
    #print(tag_to_dist.values())
    #exit()
    new_x = np.linspace(0,100,300)
    open_spline = interp1d(range(0,101),opens,kind="cubic")
    close_spline = interp1d(range(0,101),closes,kind="cubic")
    if(tag=="6AA_f12"):
      ax.plot(new_x,open_spline(new_x),'-',color='black',linewidth=1,label="original pore")
      ax.plot(new_x,close_spline(new_x),'-',color='black',linewidth=3, label="new pore")
      ax.legend(ncol=2, loc = "lower right")
    else:
      ax.plot(new_x,open_spline(new_x),'-',color='black',linewidth=1)
      ax.plot(new_x,close_spline(new_x),'-',color='black',linewidth=3)
 
  ## End of loop over tags
  #vert_shift+=6
  
#plt.legend(loc="lower right", fontsize=18)
plt.show()
fig.savefig('{}_LARGE.png'.format("ALL_Tiled"), dpi=300)
  
'''
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

'''
