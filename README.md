# DCH_detector

## Before starting

This detector uses a shape called twisted tube, which has a bug in current release of Geant4. I think I fixed it, and the patch was merged into Geant4, but I am not sure how long it will take to percolate to key4hep and other software stacks. One way to reproduce the bug is doing the following:

```
cd /tmp/
git clone -b ttbug2 git@github.com:atolosadelgado/G4_simple_box.git
cd G4_simple_box/
ls
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build -j 4 
./build/main 
```

If the output looks like the following, it means the patch was not applied yet:
```
First call: 9.2388
Second call: 1.79769e+308
Third call: 1.79769e+308
32.2836
9e+99
```

The expected behavior would lead to the following output:
```
First call: 9.2388
Second call: 9.2388
Third call: 9.2388
32.2836
32.2836
```


## To compile

Source nightlies (or run with local installation to avoid problems with the software stack)

```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

Build this project with `cmake`
```bash
cmake -B build -S . -D CMAKE_INSTALL_PREFIX=install
cmake --build build -- install
k4_local_repo
```

## To run simulation

```
ddsim --steeringFile steering.py --outputFile 'dch_proton_10GeV.root' -N 10 --runType batch
```

## To run digitization

```
wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root
k4run DCHdigi/test/runDCHdigi.py
```
