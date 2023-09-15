# HZG BDT ReadMe

# Procedure for Ntuples:
1. To make ntuples for BDTS you first need to ssh onto cms37
2. There are to run the X_ntuple_producer.py file
  a. First the correct version of root can be set with:

  source /data1/jbkim/Linux/el7_v1/root-6.24.02/bin/thisroot.sh
  
  export LD_LIBRARY_PATH=/data1/jbkim/Linux/el7_v1/lib:$LD_LIBRARY_PATH
  
  export LD_LIBRARY_PATH=/data1/jbkim/Linux/el7_v1/lib64:$LD_LIBRARY_PATH

  b. The sceond component is the specific version of python which is found here:

  /data1/jbkim/Linux/el7_v1/bin/python3.9
  
4. The ntuples can then be made with the ggF_ntuple_producer.py/VBF_ntuple_producer.py


# Procedure for making BDTs:
1. Make ntuples using procedure above
2. run the BDT producer file with root (e.g. root ggF_BDT_producer.C)

The current set-up has no input requirements, but development can happen to improve the BDT making process.
