import numpy as np
import math

#---------------paramters-------------------     
nmolmax = 1124
nmolet = 1000
natomet = 4
natommax = 4
nmolwa = 124
natomwa = 3
nstep = 8000
nreq = 50
nbin = 50
nbinmax = 30
boxl = 47.716003984042
hboxl = boxl * 0.5
zero = 0.0
delr = 0.1
const = 0.5
dlx = 0.5

posx = np.zeros((nstep, nmolmax, natommax))
posy = np.zeros((nstep, nmolmax, natommax))
posz = np.zeros((nstep, nmolmax, natommax))
midx = np.zeros((nstep, nmolmax, natommax))
midy = np.zeros((nstep, nmolmax, natommax))
midz = np.zeros((nstep, nmolmax, natommax))
m1 = np.zeros((nstep, nmolmax, natommax))

poswax = np.zeros((nstep, nmolmax, natommax))
posway = np.zeros((nstep, nmolmax, natommax))
poswaz = np.zeros((nstep, nmolmax, natommax))
m1w = np.zeros((nstep, nmolmax, natommax))

cx = np.zeros((nstep, nmolmax))
cy = np.zeros((nstep, nmolmax))
cz = np.zeros((nstep, nmolmax))
cmx = np.zeros((nstep, nmolmax))
cmy = np.zeros((nstep, nmolmax))
cmz = np.zeros((nstep, nmolmax))
comx = np.zeros((nstep, nmolmax))
comy = np.zeros((nstep, nmolmax))
comz = np.zeros((nstep, nmolmax))

cosvalue = np.zeros((nreq, nmolmax))
valtot = np.zeros((nreq, nmolmax), dtype=float)
c2tofp = np.zeros((nreq, nmolmax))

icount = np.zeros((nreq, nmolmax), dtype=int)

counterbin = np.zeros(nbinmax, dtype=int)

t = np.zeros((nbin, nmolmax), dtype=float)

finalvalue = np.zeros(nmolmax, dtype=float)
consta = np.zeros(nmolmax, dtype=float)
constb = np.zeros(nmolmax, dtype=float)
constc = np.zeros(nmolmax, dtype=float)
lt = np.zeros(nmolmax, dtype=float)

valzerotot = 0.0
conb = conc = 0.0


# --- Read HISTORY ------------------------------------------
with open("HISTORY", "r") as f:
     f.readline()
     f.readline()

     for it in range(nstep):
         f.readline()
         f.readline()
         f.readline()
         f.readline()

         for i in range(nmolmax):
             if i < nmolet:
                for j in range(natomet):
                    line = f.readline().split()
                    if len(line) < 4:
                        continue
                    m, g, mass, m2 = str(line[0]), int(line[1]), float(line[2]), float(line[3])
                    line=f.readline().split()
                    if len(line) < 3:
                        continue
                    px, py, pz = float(line[0]), float(line[1]), float(line[2])
                    posx[it, i, j], posy[it, i, j], posz[it, i, j] = px, py, pz
                    xtarget, ytarget, ztarget = posx[it, i, 0], posy[it, i, 0], posz[it, i, 0]

                    if (xtarget - px) >= hboxl: px += boxl
                    if (xtarget - px) < -hboxl: px -= boxl
                    if (ytarget - py) >= hboxl: py += boxl
                    if (ytarget - py) < -hboxl: py -= boxl
                    if (ztarget - pz) >= hboxl: pz += boxl
                    if (ztarget - pz) < -hboxl: pz -= boxl

                    f.readline()
                    f.readline()

             else:
                 for j in range(natomwa):
                     line = f.readline().split()
                     if len(line) < 4:
                        continue
                     mw, gw, massw, m2w = str(line[0]), int(line[1]), float(line[2]), float(line[3])
                     line=f.readline().split()
                     if len(line) < 3:
                         continue
                     pxwa, pywa, pzwa =float(line[0]), float(line[1]), float(line[2])
                     poswax[it, i, j], posway[it, i, j], poswaz[it, i, j] = pxwa, pywa, pzwa
                     f.readline()
                     f.readline()
#-------------------------------------------------------------------------
##with open("check1.dat","a") as fo1:
##     for it in range(nstep):
##         for i in range(nmolet):
##             for j in range(natomet):
##                fo1.write(f"{it:d} {i:d} {j:d} {posx[it, i, j]:.2f} {posy[it, i, j]:.2f} {posz[it, i, j]:.2f}\n")

#----------------for zero time gap----------------------------------------
#     for j in range(nstep):
#         for k in range(nmolet):
#             a1 = posx[j, k, 3] - posx[j, k, 2]
#             a2 = posy[j, k, 3] - posy[j, k, 2]
#             a3 = posz[j, k, 3] - posz[j, k, 2]
#             with open("try1.dat", "a") as fw:
#                      fw.write(f"{j:.2f} {k:.2f} {posx[j,k,3]:.2f} {posy[j,k,3]:.2f} {posz[j,k,3]:.2f}\n")
#
#             if (a1 >= hboxl):
#                 a1 -= boxl
#             if (a1 < hboxl):
#                 a1 += boxl
#             if (a2 >= hboxl):
#                 a2 -= boxl
#             if (a2 < hboxl):
#                 a2 += boxl
#             if (a3 >= hboxl):
#                 a3 -= boxl
#             if (a3 < hboxl):
#                 a3 += boxl
#
#             usq = ((a1 * a1) + (a2 * a2) + (a3 * a3))
#             magu = usq ** 0.5
#
#             b1 = posx[j, k, 3] - posx[j, k, 2]
#             b2 = posy[j, k, 3] - posy[j, k, 2]
#             b3 = posz[j, k, 3] - posz[j, k, 2]
#
#             if (b1 >= hboxl):
#                 b1 -= boxl
#             if (b1 < hboxl):
#                 b1 += boxl
#             if (b2 >= hboxl):
#                 b2 -= boxl
#             if (b2 < hboxl):
#                 b2 += boxl
#             if (b3 >= hboxl):
#                 b3 -= boxl
#             if (b3 < hboxl):
#                 b3 += boxl
#
#             vsq = ((b1 * b1) + (b2 * b2) + (b3 * b3))
#             magv = vsq ** 0.5
#
#             dotpdt = ((a1 * b1) + (a2 * b2) + (a3 * b3))
#             selfpdt = magu * magv
#
#             costheta = float(dotpdt / selfpdt)
#             cossqzero = costheta * costheta
#
#             valzero = (3.0 * cossqzero) - 1.0
#             valzerotot += valzero

#---------------------------------non zero time gap------------------

     for k in range(0,nmolet):
         for i in range(1,(nreq)):
             for j in range(0,(nstep-(i))):
                 icount[i, k] += 1
                 ij = i + j

                 a1 = posx[ij, k, 2] - posx[ij, k, 1]
                 a2 = posy[ij, k, 2] - posy[ij, k, 1]
                 a3 = posz[ij, k, 2] - posz[ij, k, 1]
  
                 if (a1 >= hboxl):
                    a1 -= boxl
                 if (a1 < (-hboxl)):
                    a1 += boxl
                 if (a2 >= hboxl):
                    a2 -= boxl
                 if (a2 < (-hboxl)):
                    a2 += boxl
                 if (a3 >= hboxl):
                    a3 -= boxl
                 if (a3 < (-hboxl)):
                    a3 += boxl

                 usq = ((a1 * a1) + (a2 * a2) + (a3 * a3))
                 magu = float(usq ** 0.5)

                 b1 = posx[j, k, 2] - posx[j, k, 1]
                 b2 = posy[j, k, 2] - posy[j, k, 1]
                 b3 = posz[j, k, 2] - posz[j, k, 1]

                 if (b1 >= hboxl):
                    b1 -= boxl
                 if (b1 < (-hboxl)):
                    b1 += boxl
                 if (b2 >= hboxl):
                     b2 -= boxl
                 if (b2 < (-hboxl)):
                     b2 += boxl
                 if (b3 >= hboxl):
                     b3 -= boxl
                 if (b3 < (-hboxl)):
                     b3 += boxl

                 vsq = ((b1 * b1) + (b2 * b2) + (b3 * b3))
                 magv = float(vsq ** 0.5)

                 dotpdt = ((a1 * b1) + (a2 * b2) + (a3 * b3))
                 selfpdt = float(magu * magv)

                 costheta =float(dotpdt / selfpdt)
                 cossq = float(costheta * costheta)
                 val = float((3.0 * cossq) - 1.0)
                 valtot[i, k] += val
##                 with open("check2.dat", "a") as fo2:
##                     fo2.write(f"{k:.2f} {i:.2f} {j:.2f} {icount[i,k]:.2f} {ij:.2f} {posz[ij,k,2]:.2f} {posz[ij,k,1]:.2f} {a3:.2f} \n")
##                     fo2.write(f"{k:.2f} {i:.2f} {j:.2f} {icount[i,k]:.2f} {ij:.2f} {b1:.2f} {b2:.2f} {b3:.2f} {vsq:.2f} {magv:.2f} \n")
##                     fo2.write(f"{k:.2f} {i:.2f} {j:.2f} {icount[i,k]:.2f} {ij:.2f} {dotpdt:.2f} {selfpdt:.2f} {costheta:.2f} {cossq:.2f} {val:.2f} \n")
##                     fo2.write(f"{k:.2f} {i:.2f} {j:.2f} {icount[i,k]:.2f} {ij:.2f} {val:.2f} {valtot[i,k]:.2f} \n")
##                     print(f"{i:.2f} {j:.2f} {ij:.2f} {k:.2f} {posx[ij,k,3]:.2f} {posy[ij,k,3]:.2f} {posz[ij,k,3]:.2f} \n") 
##                     print(f"{i:.2f} {k:.2f} {a1:.2f} {a2:.2f} {a3:.2f} {usq:.2f} {magu:.2f} \n")
##                     print(f"{i:.2f} {k:.2f} {}{dotpdt:.2f} {selfpdt:.2f} {cossq:.2f} {val:.2f} {valtot[i,k]:.2f} \n")

#-------------------------------reorientation lifetime of individual molecule--------------

     for j in range(0,nmolet):
         for i in range(1,nreq):
             valtot[i, j] = valtot[i, j] / (float(icount[i, j]))
             c2tofp[i, j] = const * valtot[i, j]
##             with open("check3.dat", "a") as fo3:
##                 fo3.write(f"{j:.2f} {i:.2f} {valtot[i,j]:.2f} {c2tofp[i,j]:.2f} \n")         
             #print(f"{i*0.5:.2f} {j:.2f} {valtot[i,j]:.2f} {float(icount[i,j]):.2f} {c2tofp[i,j]:.2f} \n")
         
#---------------------------------simpson integration for lifetime ----------------
     for j in range(0,nmolet):
         for i in range(1,(nbin-1)):
             k = ((i) % 2)
             if (k != 0):
               conb += c2tofp[i, j]
             else:
               conc += c2tofp[i, j]
         constb[j] = conb
         conb = 0.0
         constc[j] = conc
         conc = 0.0
         lt[j] =0.5*(((1.0 + c2tofp[(nbin-1), j]) + (constb[j] * 4.0) + (constc[j] * 2.0))/3.0)
##         with open("check4.dat", "a") as fo4:
##              fo4.write(f"{j:.2f} {i:.2f} {c2tofp[i,j]:.2f} {c2tofp[(nbin-1),j]:.2f} {constb[j]:.2f} {constc[j]:.2f} {lt[j]:.2f} \n")
##                 fo4.write(f"{j:.2f} {i:.2f} \n")
##         print(f"{j:.2f} {i:.2f} {c2tofp[i,j]:.2f} {c2tofp[(nbin-1),j]:.2f} {constb[j]:.2f} {constc[j]:.2f} {lt[j]:.2f}\n")

#--------------------binning of variables -----------------------------------
         binfordist = int(round((lt[j]) / delr))
         if (binfordist <= nbinmax):
             counterbin[binfordist] += 1
#---------------------------------------------------------------------------
with open("reor-two-dist-etoh-4.5wtwa-py.dat", "w") as fo5:   
      for i in range(nbinmax):
          delrt = delr * i
          ans = (float(counterbin[i]) / float(nmolet))
          fo5.write(f"{delrt:.4f} {ans:.8f}\n")




          
