This code is used to compute the efficiency of the ALICE MID RPCs starting from the AO2D objects and it uploads the CCDB object to the ALICE test ccdb

The input file is also inserted in the repo and in order to use this code you should enter the O2 environment with the command:

```
alenv enter O2/latest
```

Then enter root and exectue:

```
.x calculateEfficiency.C
```

Afterwards, once the object is uploaded to ccdb you can run a simple simulation of 100 events with MID only use the following command, inside the O2 environment:

```
o2-sim -g fwmugen -m MID -n 100
o2-sim-digitizer-workflow --condition-remap "http://ccdb-test.cern.ch:8080=MID/Calib/ChamberEfficiency"
```

where the --condition-remap is used to tell the digitizer code to fetch the object from the test instance of CCDB
