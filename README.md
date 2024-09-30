VMD post-processing script with custom functions to extract, transform and use quantum simulation data. To calculate quantum simulations related properties.

To use:
# load file and molecules
```

source analysis_vmd.tcl
mol new molecule1.psf
mol addfile molecule1.pdb waitfor -1
set n [molinfo top get id]

```
# select atoms
```

set sel1 [atomselect top "resname DIA"]
set sel2 [atomselect top "resname LIG"]

```
# functions use
```

align_fixedatoms $n top $fixed $fixed
cal_gofr $sel1 $sel2 0.01 7 2000 -1 datafiles/gofr.dat
distance $sel1 $sel2 top datafiles/distance.dat datafiles/histogram_distance.dat
avg_pos molid $sel1 datafiles/pos.dat
rmsd $sel molid1 molid2 datafiles/rmsd.dat

```
