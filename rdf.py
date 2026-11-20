import numpy as np

# ------------------ SIMULATION SETTINGS ------------------
nmolmax = 1124      # Total molecules (ethanol + water)
natommax = 4        # Max atoms per molecule
nmolet = 1000       # Number of ethanol molecules
natomet = 4         # Atoms per ethanol molecule
nmolwa = 124        # Number of water molecules
natomwa = 3         # Atoms per water molecule
nstep = 10          # Number of frames to read from HISTORY
nbinmax = 160       # Number of distance bins

# Simulation box dimensions (in Angstroms)
xlen = ylen = zlen = 47.716003984042
hxlen, hylen, hzlen = xlen / 2, ylen / 2, zlen / 2
vol = xlen ** 3

delr = 0.1          # Bin width (in Angstroms)
con = (4.0 / 3.0) * np.pi

# ------------------ ARRAYS ------------------
hist11 = np.zeros(nbinmax, dtype=int)
gr = np.zeros(nbinmax, dtype=float)

posx = np.zeros((nmolmax, natommax))
posy = np.zeros((nmolmax, natommax))
posz = np.zeros((nmolmax, natommax))
m1   = np.zeros((nmolmax, natommax))

cx = np.zeros(nmolmax)
cy = np.zeros(nmolmax)
cz = np.zeros(nmolmax)
cox = np.zeros(nmolmax)
coy = np.zeros(nmolmax)
coz = np.zeros(nmolmax)

# ------------------ READ TRAJECTORY & FIND CENTERS OF MASS ------------------
print("Reading HISTORY file...")

with open("HISTORY", "r") as f:
    # Skip first two header lines
    f.readline()
    f.readline()

    for it in range(nstep):
        # Skip 4 header lines per timestep
        for _ in range(4):
            f.readline()

        for i in range(nmolmax):
            if i < nmolet:  # Ethanol molecules
                Mt = 0.0
                cx[i] = cy[i] = cz[i] = 0.0

                for j in range(natomet):
                    line = f.readline().split()
                    mass = float(line)
                    m1[i, j] = mass

                    # Read atom positions
                    px, py, pz = map(float, f.readline().split())
                    posx[i, j], posy[i, j], posz[i, j] = px, py, pz

                    # Fix molecules split across box edges (Periodic Boundary Conditions)
                    xtarget, ytarget, ztarget = posx[i, 0], posy[i, 0], posz[i, 0]

                    if (xtarget - px) >= hxlen: px += xlen
                    if (xtarget - px) < -hxlen: px -= xlen
                    if (ytarget - py) >= hylen: py += ylen
                    if (ytarget - py) < -hylen: py -= ylen
                    if (ztarget - pz) >= hzlen: pz += zlen
                    if (ztarget - pz) < -hzlen: pz -= zlen

                    Mt += mass
                    cx[i] += mass * px
                    cy[i] += mass * py
                    cz[i] += mass * pz

                    # Skip velocity and force lines
                    f.readline()
                    f.readline()

                # Calculate Center of Mass for this ethanol molecule
                cox[i] = cx[i] / Mt
                coy[i] = cy[i] / Mt
                coz[i] = cz[i] / Mt

            else:  # Skip water molecules
                for _ in range(natomwa * 4):
                    f.readline()

        # ------------------ CALCULATE DISTANCES ------------------
        for i in range(nmolet):
            for k in range(i + 1, nmolet):
                dxx = cox[i] - cox[k]
                dyy = coy[i] - coy[k]
                dzz = coz[i] - coz[k]

                # Minimum image convention
                if dxx >= hxlen: dxx -= xlen
                if dyy >= hylen: dyy -= ylen
                if dzz >= hzlen: dzz -= zlen
                if dxx <= -hxlen: dxx += xlen
                if dyy <= -hylen: dyy += ylen
                if dzz <= -hzlen: dzz += zlen

                rlm = np.sqrt(dxx**2 + dyy**2 + dzz**2)
                nbin11 = int(round(rlm / delr)) + 1

                if nbin11 <= nbinmax:
                    hist11[nbin11 - 1] += 2

# ------------------ NORMALIZE g(r) ------------------
dbulk11 = float(nmolet) / vol
for nbin11 in range(1, nbinmax + 1):
    rlower = (nbin11 - 1) * delr
    rupper = rlower + delr
    deal11 = con * dbulk11 * ((rupper**3) - (rlower**3))
    gr[nbin11 - 1] = hist11[nbin11 - 1] / (nmolet * deal11 * nstep)

# ------------------ SAVE OUTPUT ------------------
output_filename = "gr-etoh-4.5wtwa-351K-py.dat"
with open(output_filename, "w") as fout:
    for nbin11 in range(1, nbinmax + 1):
        rlower = (nbin11 - 1) * delr
        fout.write(f"{rlower:.6f} {gr[nbin11 - 1]:.6f}\n")

print(f"Done! Results saved to {output_filename}")

