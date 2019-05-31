# Phased Array Antenna Analysis

This is a bunch of python scripts to simulate, analyze and plot the results of different simulations conducted in STK and CST Microwave Studio.

There are basically three different parts:

	- CST result analysis
	- STK result analysis
	- Python scripts for simulating the fields of phased array and patch antennas

## CST Result Analysis
The main script that analyzes the CST results is *CSTResultLoader.py*.
It takes the steering sweep simulation results from CST and parses it to plot Heatmaps. 
The *plotCSTHeatmap.py* can also be used but does not provide as much functionality (legacy script).

## STK Scripts
The STK scripts are rather self-explanatory according to their title.
Each of them has to be adjusted so that the path to the simulation results is correct.

## Antenna Toolbox
The antenna toolbox consists of 3 class definitions: *Antenna.py*, *RectangularArray.py* and *RectangularPatch.py*, 
where the two latter inherit important methods as the plotting of the 3D pattern or pattern cuts from *Antenna.py*.
All classes support the calculation of the normalized pattern function and antenna directivity. 

For example, an antenna with an arbitrary precomputed field can be instantiated by:
```Python 
patch = Antenna()
patch.loadPatternFromFile("sim/patch_8Ghz_efield_ro5880_0.254.txt")
```

Using the following, a rectangular array can be instantiated.
In this example, the previously created `patch` will be set as element pattern.
```python 
paa1 = RectangularArray(4,4,0.5,0.5,f,arraytype="uniform",sidelobe_level=30)
paa1.setElement(patch)
paa1.getArrayFactor(0,180,1,0,360,1,0,0)
paa1.getField(0,0)
paa1.calcDirectivity()
```

The *full_paa_comparison.py* script shows more functionality of this toolbox.
