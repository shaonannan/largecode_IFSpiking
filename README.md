# largecode_IFSpiking

# Conductance-based integrate-and-fire spiking neuron model for V1 
The main code (V1) is written in c++, and tested on  Ubuntu 16.04.

The topology connections from LGNs to V1 can be customized by re-writing and re-running the file 'thalam_cortex.py', and it will re-generate 'align_dlgnon/off/npy'.

Feedforward LGN inputs vary according to different visual stimuli, spatial locations and temporal properties of the visual stimuli and corresponding LGN responses can be customized by re-writing and re-running the file 'produce_lgn_spike_ONOFFLMI.py'. It will re-generate 'pythondata.mat', then run file 'datapytotxt.m' will obtain 'gl.txt' for the main code (V1).

other self-defined properties of V1 cells are in the form of '.txt' files, they can be re-generated by running 'ori_hyp.m'

## Get started with default configurations
First, generate external LGN inputs,
    run 'ori_hyp.m', generate clustorien/clusthyp/theta/phase.txt (properties of V1 cells);
    run 'thalam_cortex.py', generate connection map from LGN to V1;
    run 'produce_lgn_spike_ONOFFLMI.py', generate 'pythondata.mat'
    run 'datapytotxt.m', generate 'gl.txt'.
    
Put all the .txt files needed in your Workplace_directory, and 

    run 'g++ Pythontest_onoffLMI_1207_Thesis.cpp -o SELF_DEFINED_NAME'
    run './SELF_DEFINED_NAME'

Simulating results in Workplace_directory/... and save as .txt files.

