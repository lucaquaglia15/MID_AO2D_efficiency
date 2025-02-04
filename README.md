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


The format of the ccdb file is to have a `std::vector` of `o2::mid::ChEffCounter` objects.

Each element of this vector is a struct that contains the following information: `deId` (the ID number of the RPC within the MID detector), `columnId` (the column ...), and `lineId` (...) as well as a `std:array` of 4 elements, where element 0 is the counts on the bending plane, element 1 is the counts on the non-bending plane, element 2 is the counts on both panes, anf element 3 is the total counts.

The following figure shows the scheme of the ALICE MID with the detector numbering scheme adopted in ALICE O2 and it shows the values of `deId`

![Image](https://github.com/user-attachments/assets/5f99cb21-8148-4185-bc81-34cc96f9aebc)
