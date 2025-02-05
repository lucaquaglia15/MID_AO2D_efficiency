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

Each element of this vector is a struct that contains the following information: `deId` (the ID number of the RPC within the MID detector), `columnId` (the column of the MID i.e. the ideal vertical lines that follow the local board segmentation), and `lineId` (the line of the MID i.e. the horizontal lines that follow the local board segmentation) as well as a `std:array` of 4 elements, where element 0 is the counts on the bending plane, element 1 is the counts on the non-bending plane, element 2 is the counts on both panes, anf element 3 is the total counts.

The following figure shows the scheme of the ALICE MID with the detector numbering scheme adopted in ALICE O2 (`deId`) as well as the internal name of the detectors (MT xx IN/OUT y) where IN and OUT refer to the detectors inside or outside the LHC circumference.

![Image](https://github.com/user-attachments/assets/5f99cb21-8148-4185-bc81-34cc96f9aebc)

The following figure shows instead a scheme of the MID segmentation in the different RPCs and LBs. Notice that the columns go from 0 to 7 and they increase from the center of the MID outwards while the lines go from bottom to top. Notice that while column number increases independently of the RPCs, the lines start from 0 in each RPC (and their maximum depends on the given RPC but at the most it is 3). After the picure there are some examples to make it clear.

![Image](https://github.com/user-attachments/assets/be8f66b1-efed-4a24-bc8d-6de63175e333)

Note that any given LB reads out the space "space" region on all 4 MID planes and there is a total of 234 LB. We can know whether a muon track passed through any given LB on any plane and for this reason we number the LB from 1 to 936 (234*4) so that we also have the information on the plane.

For example, the LB identified by the number 28 is numbered 28 on MT11 while on MT12 it is numbered 28+234 = 262 and so on and so forth. 

Few examples of column/line numbering of LB:
- LB 168 on MT11 has column = 3 and line = 3 
- LB 87 on MT11 has column = 4 and line = 1
- LB 113 on MT11 has column = 7 and line = 0

This is important because each element of the `std::vector` of `o2::mid::ChEffCounter` objects must be identified by the deId, columnId and lineId values described earlier.

The CCDB objects cannot be uploaded directly by the user to the [ALICE production CCDB](http://alice-ccdb.cern.ch/browse/MID/Calib/ChamberEfficiency) but only on a local instance of the CCDB. In order to have them uploaded centrally, one has to open a [JIRA ticket](https://its.cern.ch/jira/projects/O2/summary), such as [this one](https://its.cern.ch/jira/browse/O2-5759) to have them uploaded. To follow this route, one needs to save the CCDB objects in a .root format and upload them to JIRA, together with an upload script. To this end, two codes have been developed:
