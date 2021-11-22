import numpy as np
import signal_procsssing as sp 
import f90nml
import directory_edit as de 
import phys_cons as pc 

de.make_dir("fig")

c = pc.c
params = f90nml.read("param.inp")
nx = params["space"]["nxx"]
dx = params["space"]["dx"]
dt = params["time"]["deltat"]/(c*np.sqrt(1.0/dx/dx))
nstep = params["time"]["nstep"]+1
t = np.arange(0,nstep*dt,dt)
freq = np.fft.fftfreq(nstep,dt)
freq = freq[0:int(nstep/2)]

iex = sp.input_signal("exin.txt")
iey = sp.input_signal("eyin.txt")
iez = sp.input_signal("ezin.txt")

rex = sp.input_signal("ex_ref.txt")
rey = sp.input_signal("ey_ref.txt")
rez = sp.input_signal("ez_ref.txt")

fiex = sp.fft(iex)
fiey = sp.fft(iey)
fiez = sp.fft(iez)

frex = sp.fft(rex)
frey = sp.fft(rey)
frez = sp.fft(rez)

sp.plot_signal(t,iey,color="blue",label="ey")
sp.show(draw=True,save=False,save_name="input_signal")

sp.plot_signal(freq,fiey,color="blue",label="fey",xlabel="frequency")
sp.show(draw=True,save=False,save_name="fft-input_signal")

fie = fiex**2+fiey**2+fiez**2
fre = frex**2+frey**2+frez**2

rcs = sp.rcs(fie,fre)
sp.plot(freq,rcs,color="blue",label="reflection_cof",ymin=-50,ymax=0,title="reflection-cof",xlabel="frequency")
sp.show(draw=True,save=False,save_name="reflection_cof")