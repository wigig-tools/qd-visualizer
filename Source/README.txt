This repository contains a python application (using Mayavi library) to visualize ns-3 IEEE 802.11ad with Q-D channel simulation results.

**************** Features ****************
- This application allows to visualizer the beamforming training evolution, and in particular, the SLS TxSS and RxSS phase. Each beamforming training is visualized thanks to the directivity of the antenna pattern resulting of the beamforming training SLS phase.
- The Multi-Paths components between a pair of transmitter/receiver are displayed.
- The system level performance (throughput/SNR) can be visualized

The Q-D visualizer comes with three different scenarios corresponding to the one presented in WNGW19:

- L-Shaped Room
- Spatial Sharing
- Dense deployment

**************** Usage ****************

**** L-Shaped Room scenario *****
This scenario shows the SLS results and throughput/SNR results for the communication between a STA and an AP. 
To launch the visualization:

python QDVisualizer.py -f L-Room

Then, you have to selected to set the "Sta training" to 2. To visualize the evolution, just update the "Trace index".

**** Spatial Sharing ******
This scenario shows the SLS results evolution when two BSSs communicate simultaneously. 

To launch the visualization:
python QDVisualizer.py -f SpatialSharing

Then, select the STA for whom you want to observe the SLS evolution (3 or 4) and use the "Trace index" slider to visualize its evolution. 


**** Denser Scenario  ******
This scenario shows the SLS results when 10 STAs perform their SLS beamforming training with 1 AP. 
To launch the visualization:
python QDVisualizer.py -f DenserScenario

Just select the STA you want to observe with "Sta training".





