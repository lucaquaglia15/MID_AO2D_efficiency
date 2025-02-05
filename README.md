This code is used to compute the efficiency of the ALICE MID RPCs starting from the AO2D objects and it uploads the CCDB object to the ALICE test ccdb

The input file is also inserted in the repo and in order to use this code you should enter the O2 environment with the command:

```
alienv enter O2/latest
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

The CCDB objects are produced on a run-by-run basis and they cannot be uploaded directly by the user to the [ALICE production CCDB](http://alice-ccdb.cern.ch/browse/MID/Calib/ChamberEfficiency) but only on a local instance of the CCDB. In order to have them uploaded centrally, one has to open a [JIRA ticket](https://its.cern.ch/jira/projects/O2/summary), such as [this one](https://its.cern.ch/jira/browse/O2-5759) to have them uploaded. To follow this route, one needs to save the CCDB objects in a .root format and upload them to JIRA, together with an upload script. To this end, the following files can be used:

- One needs to download the Hyperloop output file from the ALICE grid. Open the `config.yml` file and change the `output_dir_name` to your desired output, modify the `run_list` and the `alien_input_path` lists with the run numbers and the path of the files on the ALICE grid respectively. One this is done do the following to download the data from the ALICE grid, creating one folder for each run and each folder contains the `AnalysisResults.root` file inside.

```
$ alienv enter O2Physics/latest
$ python utils.py config.yml —download
```

- Once the files have been downloaded, one needs to fetch the information on each run (i.e. start/end time, interaction rate, magnetic field polarity etc..). To do this, you can use the `getInteractionRate.py` script (that has to be in the same folder as `task.C`). Open it and modify the list `runs` to include the same runs you have downloaded from the ALICE grid and modify also the `basePath`/`period` variables to point in the same folder as the one containg all the runs you just downloaded. Once done, launch the script with the following command. This will create a .txt file called `run_IR_Bfield.txt` that has different columns (in order: 0/1 if the run has to be considered or not >> run number >> IR (it makes sense only in Pb-Pb) >> megntic field polarity >> start of run >> end of run. This file is then processed further by other scripts.

```
$ python3 getInteractionRate.py python3 getInteractionRate.py -f 551418 -l 551427 --run mdquality --duration 1
```
- The first one is `produceObjects.C` and the second one is `prepareCCDBUpload.C`. Before running them you should modify the paths in both file to point to the same folder as where all the files from all the runs have been saved. Then do the following (this will create one .root file for each run and it will also create the upload sciprt `upload.sh` in the same folder as where all the runs have been saved).
```
$ alienv enter O2Physics/latest
$ root produceObjects.C
$ root prepareCCDBUpload.C
```
- When this is finished you just have to zip the folder where the .root files are saved and upload it to JIRA.

Here are some notes on how the CCDB objects are saved to a .root file. This is done using the `WriteObjectAny` member function of the `TDirectory` ROOT CERN class which takes as first argument the pointer to the object (in this case the `std::vector` of `o2::mid::ChEffCounter` objects) and the second argument is the object type (`std::vector<o2::mid::ChEffCounter>`) and the third element is the name with which they will be saved in the .root file (`ccdb-object`).

Once the obejects are saved in a .root file, one can check the content in this way (this are lines of code to be launched inside of root in the O2Physics environment but they can be included in a macro).

```
$ alienv enter O2Physics/latest
$ root
$ #
$ TFile *f = new TFile("o2-mid-ChEffCounter_xxxxxx.root","READ") //xxxxxx is the run number
$ f->cd()
$ vector<o2::mid::ChEffCounter> *v = (vector<o2::mid::ChEffCounter>*)f->Get("ccdb-object")
$  cout << v->at(yy).getCounts(o2::mid::EffCountType::BothPlanes) //to print out the value of the specific counter (yyy is the vector element
```
- The vector element depends on how the vector was filled. In this case I have used three loops. A first one on the detector elements (0-71), then one on the columns (0-7) and one on the lines (0-depends on the specific deId). Referring to the LB segmentation scheme (to give an idea), the elements are as follows: 
	- v->at(0) corresponds to LB1
	- v->at(1) corresponds to LB17
	- v->at(2) corresponds to LB39
	- v->at(3) corresponds to LB61
	- v->at(4) corresponds to LB77 
	- v->at(5) corresponds to LB93
	- v->at(6) corresponds to LB109
	
	- v->at(15) corresponds to LB78
	
	- v->at(34) corresponds to LB7
	
and so on an so forth 


