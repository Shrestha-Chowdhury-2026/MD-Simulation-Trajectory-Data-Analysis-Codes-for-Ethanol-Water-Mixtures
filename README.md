# Ethanol Center-of-Mass RDF Calculator

This Python script reads a DL_POLY `HISTORY` trajectory file of an ethanol–water mixture and calculates the center-of-mass Radial Distribution Function, $g(r)$, for ethanol molecules.

## How to Run
1. Make sure you have Python and `numpy` installed (`pip install numpy`).
2. Place your DL_POLY `HISTORY` file in the same folder as `rdf.py`.
3. Run the script:
   ```bash
   python rdf.py
   ```
4. The output will be saved as `gr-etoh-4.5wtwa-351K-py.dat`.

# Ethanol Center-of-Mass MSD Calculator

This Python script reads a DL_POLY `HISTORY` trajectory file of an ethanol–water mixture and calculates the center-of-mass Mean Square Displacement, $MSD$, for ethanol molecules.

## How to Run
1. Make sure you have Python and `numpy` installed (`pip install numpy`).
2. Place your DL_POLY `HISTORY` file in the same folder as `MSD.py`.
3. Run the script:
   ```bash
   python MSD.py
   ```
4. The output will be saved as `MSD-etoh-4.5wtwa-py.dat`.

