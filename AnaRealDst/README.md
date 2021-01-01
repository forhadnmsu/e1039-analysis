# e1039-analysis/AnaRealDst

A compact program to analyze the "DST" file, which is created by the Main DAQ decoder.
Since it is yet under development and also depends on other software under development,
it might stop working occasionally.  In such cases, please contact the author.

## Usage:

The primary uses of this module is to do high voltage optimizations for all the hodoscopes. Using (x,y,z) hit positions from different detector comsic muon/single tracks have been created. Using these tracks efficiency of individual paddles can be studied. To make sure we do not have any biased in the analysis study, we also have applied the hit conditions, so efficiency-fining plane will not be involved in the track building, and even there was any contribtuions, these tracks won't be taken acount to find the efficiency numbers.  

Since different run number represents different conditions (especially the high voltages), we have tracked the run numbers and the corresponding HV of hodoscopes, so we can easily se the plateau curve.


1. source setup.sh
1. cd work
1. cmake-e1039-ana
1. make install
1. root -b -q '../macro/Fun4RealDst.C(202)'

Here are some tips:

- You need execute "source setup.sh" every time you open a new terminal.
- You need execute "make install" every time you modify any file in "src/".
- You can execute the macro at any directory.

## Advanced Usage

- The contents of the analysis are defined in class "AnaRealDst".
  You can modify this class (i.e. "src/AnaRealDst.cc") to implement your own analysis.
    - The E1039 data are structured by the SQ interface classes.  To access them, you use several member functions of these classes.  You can refer to the [oline document](https://e1039-collaboration.github.io/e1039-doc/annotated.html).  Below are frequently-used classes;
        - [SQRun](https://e1039-collaboration.github.io/e1039-doc/d7/db7/classSQRun.html)
        - [SQEvent](https://e1039-collaboration.github.io/e1039-doc/d9/dd7/classSQEvent.html)
        - [SQHitVector](https://e1039-collaboration.github.io/e1039-doc/d9/dbc/classSQHitVector.html)
        - [SQHit](https://e1039-collaboration.github.io/e1039-doc/de/d79/classSQHit.html)
        - [GeomSvc](https://e1039-collaboration.github.io/e1039-doc/d0/da0/classGeomSvc.html)
- The run ID and the number of events for analysis are defined in "macro/Fun4RealDst.C".
  You can change them so that you need not give them in the command line.

## Author

Kenichi Naknao <knakano@nucl.phys.titech.ac.jp>
