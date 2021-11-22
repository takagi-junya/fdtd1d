#------------------------------------------
# Plot 3D data in 1D
#------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import h5py
import os
#import gc
import directory_edit as de
import sys
import matplotlib.patches as patches
import shutil
import make_mp4
import f90nml
import tkinter
from tkinter import ttk 
  
  
def init():
    params = f90nml.read("param.inp")
    nx = params["space"]["nxx"]
    amp = np.max(params["wave"]["amps"])
    xmin = 0 
    xmax = nx 
    ymin = -0.25*amp 
    ymax = 0.25*amp
    return xmin,xmax,ymin,ymax

def plot1D(data,vmin,vmax,xmin,xmax,flag,area):
  if(flag==0):
    img = plt.plot(data,color="red")
  else:
    img = plt.plot(data,color="red")
  if(area==1):
    plt.xlim([xmin,xmax])
    
  return img

def vis(comp,type,area,interval,xmin,xmax,ymin,ymax):
  k = 0
  keys=[]
  ims = []
  cwd = os.getcwd()
  os.chdir(cwd)
  print("comp:",comp)
  print("type:",type)
  print("area:",area)
  
  z0   =  376.73031
  fig,ax = plt.subplots()
  
  if(type=="gif"):
    if not de.find("gif"):
        de.make_dir("gif")
      
  if(type=="png" or type=="mp4"):
    name = "1dfig_"+comp
    filename = "1dfig_"+comp
    if not de.find(name):
        de.make_dir(name)
      
  with h5py.File(comp+".h5",mode="r") as f:
    for i in f[comp].keys():
      keys.append(i)

    for k in range(0,len(keys)):
      if(k%interval!=0):
        continue;
      print(keys[k])
      lines = []
      if(type=="png" or type=="mp4"):
        fig,ax = plt.subplots()

      #data = np.array(f[comp][str(keys[k])][py,:])
      data = np.array(f[comp][str(keys[k])])
      if(comp in ["ex","ey","ez"]):
        if(area=="normal"):
          im = plot1D(data=data,vmin=ymin,vmax=ymax,xmin=xmin,xmax=xmax,flag=0,area=0)
        else:
          im = plot1D(data=data,vmin=ymin*0.5,vmax=ymax*0.5,xmin=xmin,xmax=xmax,flag=0,area=1)
      elif(comp in ["hx","hy","hz"]):
        if(area=="normal"):
          im = plot1D(data=data,vmin=ymin/z0,vmax=ymax/z0,xmin=xmin,xmax=xmax,flag=0,area=0)
        else:
          im = plot1D(data=data,vmin=ymin*0.5/z0,vmax=ymax*0.5/z0,xmin=xmin,xmax=xmax,flag=0,area=1)
      else:
        if(area=="normal"):
          im = plot1D(data=data,vmin=ymin*0.01,vmax=ymax*0.01,xmin=xmin,xmax=xmax,flag=0,area=0)
        else:
          im = plot1D(data=data,vmin=ymin*0.01,vmax=ymax*0.01,flag=0,area=1)

      plt.title(comp+" on y=1280  "+str(keys[k]).zfill(4))
      lines.extend(im)
      ims.append(lines)
      plt.grid()
      if(type=="png" or type=="mp4"):
        plt.savefig(filename+"/"+str(keys[k]).zfill(4)+".png")
        plt.close()
    if(type=="gif"):
      os.chdir("./gif/")
      ani = animation.ArtistAnimation(fig,ims,interval=250,blit=False)
      ani.save(comp+"1dx.gif",writer="imagemagick")
    if(type=="mp4"):
      print(filename)
      os.chdir(cwd)
      os.chdir(filename)
      make_mp4.mp4(comp+"_"+area)

args = sys.argv[1:]
comp = "ey"
type =  "gif"
area =  "normal"
interval = 10
xmin,xmax,ymin,ymax = init()

vis(comp,type,area,interval,xmin,xmax,ymin,ymax)