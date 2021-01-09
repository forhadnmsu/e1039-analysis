# e1039-analysis/AnaRealDst (HV-Efficiency Branch)
The primary uses of this module is to do high voltage optimizations for all the hodoscopes. This one has been developed from the basic detector analysis module [AnaRealDst](https://github.com/E1039-Collaboration/e1039-analysis), and has been extended and further developed for the high-voltage-scanning and the efficiency study.

## Produces to get histograms and ntuple.

1. source setup.sh
1. cd work
1. cmake-e1039-ana
1. make install
1. root -b -q '../macro/Fun4RealDst.C(202)'

use Fun4MultiRealDst.C for multiple runs.

Here are some tips:

- You need execute "source setup.sh" every time you open a new terminal.
- You need execute "make install" every time you modify any file in "src/".
- You can execute the macro at any directory.

## Track Building and Getting the efficiency plots

1. Using (x,y,z) hit positions from different detectoris, comsic muon/single tracks can be formed. Using these tracks, efficiency of individual paddles can be studied. Add all the detector planes (x, and y) for the track building using [UtilHodo2::Track2D](https://github.com/forhadnmsu/e1039-analysis/blob/hodo-hvScan/AnaRealDst/src/UtilHodo2.h).
1. You have to make sure that the track building detectors alwas have a hit. 
1. Since, you will have a staright line fitting, the lower chi^2 will give you better track quality, but you may add additional track quality requirement to make sure that the effciency number is better, but at the same time one need to be careful so we do not introduces any biases for the efficiency plots.

## HV- Scanning
For High voltage scanning, it is very liklely you will use a lot of different DAQ runs, and these runs could have different hodoscope high volatges. It's entrire up to you how you want to coordinate the high volages, with the DAQ run. For my study, I have took notes the high volatges with the different DAQ run numbers and associated them in the histogram paddles. You can do this way or you follow the way you prefer.  

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
Regarding the Efficiency Analysis: Forhad Hossain <forhad16@nmsu.edu>
Regarding Core Software Developement: Kenichi Naknao <knakano@nucl.phys.titech.ac.jp>
