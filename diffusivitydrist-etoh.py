import numpy as np
#---------------paramters-------------------     
nmolmax = 1124
nmolet = 1000
natomet = 4
natommax = 4
nmolwa = 124
natomwa = 3
nstep = 2000
nreq = 300
gap = 100
nbinmax = 100
starttime = 50
boxl = 47.716003984042
hboxl = boxl * 0.5
zero = 0.0
six = 6.0
timestepvalue = 0.5
delr = 0.085

posx = np.zeros((nstep, nmolmax, natommax))
posy = np.zeros((nstep, nmolmax, natommax))
posz = np.zeros((nstep, nmolmax, natommax))
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

xpre = np.zeros((nmolmax))
ypre = np.zeros((nmolmax))
zpre = np.zeros((nmolmax))

usq = np.zeros((nmolmax, nreq), dtype=float)

icount = np.zeros((nmolmax, nreq), dtype=int)
counterbin = np.zeros(nbinmax, dtype=int)

# --- Read HISTORY ------------------------------------------
with open("HISTORY", "r") as fin:
     fin.readline()
     fin.readline()

     for it in range(nstep):
         fin.readline()
         fin.readline()
         fin.readline()
         fin.readline()

         for i in range(nmolmax):
             if i < nmolet:
                Mt = 0.0
                cx[it, i] = 0.0
                cy[it, i] = 0.0
                cz[it, i] = 0.0
                xtarget = ytarget = ztarget = 0.0

                for j in range(natomet):
                    line = fin.readline().split()
                    if len(line) < 4:
                        continue
                    m, g, mass, m2 = str(line[0]), int(line[1]), float(line[2]), float(line[3])
                    m1[it, i, j] = mass
                    line=fin.readline().split()
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

                    Mt += m1[it, i, j]
                    cx[it, i] += mass * px
                    cy[it, i] += mass * py
                    cz[it, i] += mass * pz

                    fin.readline()
                    fin.readline()

                comx[it, i] = cx[it, i] / Mt
                comy[it, i] = cy[it, i] / Mt
                comz[it, i] = cz[it, i] / Mt
             else:
                 Mt = 0.0
                 for j in range(natomwa):
                     line = fin.readline().split()
                     if len(line) < 4:
                        continue
                     mw, gw, massw, m2w = str(line[0]), int(line[1]), float(line[2]), float(line[3])
                     m1w[it, i, j] = massw
                     line=fin.readline().split()
                     if len(line) < 3:
                         continue
                     pxwa, pywa, pzwa =float(line[0]), float(line[1]), float(line[2])
                     poswax[it, i, j], posway[it, i, j], poswaz[it, i, j] = pxwa, pywa, pzwa
                     fin.readline()
                     fin.readline()


#---------------------------unfolding of com-------------------------

     for i in range(nmolet):
        xsh = 0.0
        ysh = 0.0
        zsh = 0.0
        for it in range(nstep):
            cmx[it, i] = 0.0
            cmy[it, i] = 0.0
            cmz[it, i] = 0.0
            if it != 0:
                xd = yd = zd = 0.0
                if (abs((comx[it, i]) - (xpre[i]))) >= hboxl:
                    xd = (comx[it, i]) - (xpre[i])
                if xd > 0.0:
                   xsh -= boxl
                if xd < 0.0:
                   xsh += boxl

                if (abs((comy[it, i]) - (ypre[i]))) >= hboxl:
                    yd = (comy[it, i]) - (ypre[i])
                if yd > 0.0:
                   ysh -= boxl
                if yd < 0.0:
                   ysh += boxl

                if (abs((comz[it, i]) - (zpre[i]))) >= hboxl:
                    zd = (comz[it, i]) - (zpre[i])
                if zd > 0.0:
                   zsh -= boxl
                if zd < 0.0:
                   zsh += boxl

            cmx[it, i] = comx[it, i] + xsh
            cmy[it, i] = comy[it, i] + ysh
            cmz[it, i] = comz[it, i] + zsh

            xpre[i] = comx[it, i]
            ypre[i] = comy[it, i]
            zpre[i] = comz[it, i]
# --- MSD calculation -------------------------------------------------
     for k in range(nmolet):
         for i in range(0, nreq):
             for j in range (0, (nstep - (i))):
                  icount[k, i] += 1
                  ij = i + j
                  
                  s1 = (((cmx[ij, k]) - (cmx[j, k])) ** 2)
                  s2 = (((cmy[ij, k]) - (cmy[j, k])) ** 2)
                  s3 = (((cmz[ij, k]) - (cmz[j, k])) ** 2)

                  u2 = ((s1) + (s2) + (s3))
                  usq[k, i] += u2
#                  with open("ufold-py.dat","a") as fo4:
#                      fo4.write(f"{ij:d} {j:d} {k:d} {u2:.2f} {usq[i]:.2f}\n")
#------------------------------------------------------------------------
     for k in range(nmolet):
         for i in range (nreq):
             usq[k, i] = usq[k, i]/(float(icount[k, i]))
#----------------------------------------------------------------------
     diffusivity = 0.0
     totalavg = 0.0
     for k in range(nmolet):
          diffusivity = 0.0
          ii = 0

          for i in range(starttime,nreq):
              if(((i-starttime) % gap)) == 0.0:
                  if (i > starttime):
                      diff = usq[k, i] - usq[k, (i-gap)]
                      slope = diff/(gap*timestepvalue)
                      diffusivity += (slope/(float(6)))
                      ii += 1
          avgdiffusivity = diffusivity/float(ii)
          binval = int(round(avgdiffusivity/delr))
          if (binval <= nbinmax):
             counterbin[binval] += 1
          totalavg += avgdiffusivity

     print(f"{(totalavg/float(nmolet)):.4f}\n")
#--------------------------------------------------------------
with open("D-dist-etoh-4.5wtwa-py.dat", "w") as fw:
     for i in range(0,nbinmax):
         t = (float(i)) * delr
         ans = (float(counterbin[i]) / float(nmolet))
         fw.write(f"{t:.2f} {ans:.16f}\n")

#-------------------------------------------------------------
     v = 0.0
     for i in range(0,nbinmax):
         v += ((float(i)) * delr) * (float(counterbin[i])/(float(nmolet)))
     print(f"{v:.4f}\n")

