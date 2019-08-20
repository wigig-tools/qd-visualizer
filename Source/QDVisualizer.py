# First, and before importing any Enthought packages, set the ETS_TOOLKIT
# environment variable to qt4, to tell Traits that we will use Qt.
import os
import pyqtgraph as pg
os.environ['ETS_TOOLKIT'] = 'qt4'
# By default, the PySide binding will be used. If you want the PyQt bindings
# to be used, you need to set the QT_API environment variable to 'pyqt'
#os.environ['QT_API'] = 'pyqt'

# To be able to use PySide or PyQt4 and not run in conflicts with traits,
# we need to import QtGui and QtCore from pyface.qt
from pyface.qt import QtGui, QtCore
from PyQt4.QtGui import *
from PyQt4.QtCore import *
# Alternatively, you can bypass this line, but you need to make sure that
# the following lines are executed before the import of PyQT:
#   import sip
#   sip.setapi('QString', 2)

from traits.api import HasTraits, Instance, on_trait_change, Range, \
   Button, Bool, Enum, List
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor
import os
import cmath 
from mayavi import mlab
from mayavi.api import Engine
from mayavi.mlab import *
import csv
import math
from mayavi.modules.labels import Labels
from mayavi.tools.mlab_scene_model import \
    MlabSceneModel
from traitsui.api import View, Item, HGroup, VGroup, Group, HSplit
from numpy import linspace, pi, cos, sin
from mayavi.core.ui.mayavi_scene import MayaviScene
import vtk
from tvtk.api import tvtk
import numpy as np
import argparse
from pyface.api import GUI

#####################################
#### Q-D Visualization software #####
#####################################
# Author: Tanguy Ropitault          #
# email: tanguy.ropitault@nist.gov  #
#####################################

######################################################################################################
# NIST-developed software is provided by NIST as a public service. You may use, copy                 #
# and distribute copies of the software in any medium, provided that you keep intact this            #
# entire notice. You may improve, modify and create derivative works of the software or              #
# any portion of the software, and you may copy and distribute such modifications or                 #
# works. Modified works should carry a notice stating that you changed the software                  #
# and should note the date and nature of any such change. Please explicitly                          #
# acknowledge the National Institute of Standards and Technology as the source of the                #
# software.                                                                                          #
#
# NIST-developed software is expressly provided "AS IS." NIST MAKES NO                               #               
# WARRANTY OF ANY KIND, EXPRESS, IMPLIED, IN FACT OR ARISING BY                                      #
# OPERATION OF LAW, INCLUDING, WITHOUT LIMITATION, THE IMPLIED                                       #
# WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,                                     #
# NON-INFRINGEMENT AND DATA ACCURACY. NIST NEITHER REPRESENTS                                        #
# NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE                                            #
# UNINTERRUPTED OR ERROR-FREE, OR THAT ANY DEFECTS WILL BE                                           #
# CORRECTED. NIST DOES NOT WARRANT OR MAKE ANY REPRESENTATIONS                                       #
# REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,                                          #
# INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY,                                            #
# RELIABILITY, OR USEFULNESS OF THE SOFTWARE.                                                        #
#                                                                                                    #
#                                                                                                    #
# You are solely responsible for determining the appropriateness of using                            #
# and distributing the software and you assume all risks associated with its use, including          #
# but not limited to the risks and costs of program errors, compliance with applicable               #
# laws, damage to or loss of data, programs or equipment, and the unavailability or                  #
# interruption of operation. This software is not intended to be used in any situation               #
# where a failure could cause risk of injury or damage to property. The software                     #
# developed by NIST employees is not subject to copyright protection within the United               #
# States.                                                                                            #
######################################################################################################


######################################################################################################
###############                          GLOBAL VARIABLE                               ###############
######################################################################################################

colorMapCustom = 'blue-red'
#########################################
######    Scenario parameters       #####
#########################################
# These parameters are obtained from configFile file
numberOfAp = 1
numberOfSta = 1
nbSectorsPerAntenna = 37 
nbAntennaAp = 1 
nbAntennaSta = 1 
azimuthCardinality = 361
elevationCardinality = 181
graphData = 0

#########################################
######    Visualizer parameters     #####
#########################################
desiredReflectionOrder = 1 # Display reflection order up to this value
trace = "" # The trace to display (based on the Mobility generated by the Q-D Realization software)

# Folder variables
scenarioFolder = "" # The folder containing all the files needed to display a given scenario
positionFile = "NodePositionsTrc" # Prefix for nodes positions (STAs + AP)
positionFolder = "NodePositions/" # The folder containing the nodes positions
mpcFolder = "MpcCoordinates/" # Folder that contains the MPC coordinates files
resultsFolder = "Results/" # Folder containing the results files
CodebookFolder = "../Codebook/" # Folder containing the different codebook files
topologyFolder = "RoomCoordinates/" # The folder containing the file defining the environment geometry

# File variable
configFile = "config.csv" # The file containing the parameters of the scenario
slsResultsFile = "slsResults.csv" # The file containing the SLS phase results 
antennaArrangementFile = "antennaPosition28.csv" # File containing the coordinates of the antenna elements
CodebookFileName = "CODEBOOK_URA_AP_28x.txt" # The codebook file used for the simulations in ns-3
topologyFilename = "RoomCoordinates.csv" # The environment geometry file
throughputMCS9File = 'mcs9.csv' # The throughput file for MCS9 (L-ROOM scenario only)
throughputMCS12File = 'mcs12.csv' # The throughput file for MCS12 (L-ROOM scenario only)
snrFile  = "snr.csv" # The SNR file f(L-ROOM scenario only)

# The data needed to represent the antenna pattern for STAs
PATTERN_STA = []
# The data needed to represent the antenna pattern for AP
PATTERN_AP = []

#########################################
##########       View 1       ###########
#########################################
# View 1 Global variable
STA_BEST_PATTERN_VIEW1 = [] # A list containing the STA 3D response patterns to display on the View 1
AP_BEST_PATTERN_VIEW1 = [] # A list containing the AP 3D response patterns to display on the View 1 
AP_VIEW1 = 0 # Used to display the APs nodes
STA_VIEW1 = 0 # Used to display the STAs nodes

#########################################
##########       View 2       ###########
#########################################
ANTENNA_PATTERN_AP_VIEW2 = [] # A list containing the AP 3D response patterns to display on the View 2
ANTENNA_PATTERN_STA_VIEW2 = [] # A list containing the STA 3D response patterns to display on the View 2
MPC = [] # A list containing the MPC to display on the View 2
ANTENNA_GEOMETRY_AP_VIEW1 = [] # APs antenna array element arrangement 
ANTENNA_GEOMETRY_STA_VIEW1 = [] # STAs antenna array element arrangement 

######################################################################################################
##############                      FUNCTIONS TO LOAD INPUT FILES                      ###############
######################################################################################################

# Call all functions related to files parsing and handle command-line parsing
def loadInput():
    global scenarioFolder
    global quality
    global apAntennaGeometry
    global staAntennaGeometry
    global nbTraces
    global graphData
    #####################################################
    ###########     Command-line parsing         ########
    #####################################################

    parser = argparse.ArgumentParser()
    # By default, we force the user to provide a scenario folder
    parser.add_argument('-f', action='store', dest='scenarioFolder',
                        help='The scenario folder')
    # By default, the first trace is displayed, i.e, beginning of the simulation
    parser.add_argument('-t', nargs = '?', action='store', dest='trace',
                        help='The Trace to visualize', default = "0")
    # By default, we order or reflection to display to 1
    parser.add_argument('-r', nargs = '?', action='store', dest='reflectionOrder',
                        help='The Desired Order of Reflection to display for the MPC',const = 1 , type = int, default = 1)
    # By default, we set the antenna pattern quality to the maximum
    parser.add_argument('-q', nargs = '?', action='store', dest='qualityOfPattern',
                        help='The Desired Quality when displaying the antenna pattern (1=Max 10=Worst)', default = 5 , type = int)
    argument = parser.parse_args()
    scenarioFolder = argument.scenarioFolder + "/"
    trace = argument.trace
    desiredReflectionOrder = argument.reflectionOrder
    quality = argument.qualityOfPattern

    #####################################################
    ##########        Load input files           ########
    #####################################################
    readConfigFile() # Read the scenario configuration
    apAntennaGeometry = []
    for i in range(numberOfAp):
        apAntennaGeometry.append(readAntennaPosition()) # Read antenna element position for the APs

    staAntennaGeometry = []
    for i in range(numberOfSta):
        staAntennaGeometry.append(readAntennaPosition()) # Read antenna element position for the STAs
   
    apCoordinates, staCoordinates, nbTraces = readAllPositionFile() # Read all nodes STAs + APs) position throughout the simulation
    if graphData:
        readThroughputResults(True) # Read throughput values
        readThroughputResults(False) # Read throughput values
        readSnrResults()
    readAllSLSResults() # Read SLS phase results
    loadCodebook()


#########################################
########    SCENARIO CONFIG  ############
#########################################
# Load the configuration of the scenario (This file is manually created)
def readConfigFile():
    global numberOfAp, numberOfSta, nbSectorsPerAntenna , nbAntennaAp,  nbAntennaSta, graphData
    global scenarioFolder
    filename = scenarioFolder +  configFile
    try: 
        with open(filename) as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                numberOfAp = int(row['NB_AP'])
                numberOfSta =   int(row['NB_STA'])
                nbSectorsPerAntenna =   int(row['NB_SECTOR'])
                nbAntennaAp =   int(row['NB_ANTENNA_AP'])
                nbAntennaSta =   int(row['NB_ANTENNA_STA'])
                graphData = int(row['GRAPH_DATA'])
    except FileNotFoundError:
        print("Error: Scenario config file: ", filename, " does not exist")
        exit

#########################################
########## MPC COORDINATES    ###########
#########################################
# Read MPC file coordinate for a given receiver/transmitter pair (Provided by Q-D realization software)
def readMPCCoordinateFile(transmitter, receiver, reflectionOrder,traceIndex):
    filename = scenarioFolder + mpcFolder + "MpcTx" + str(transmitter) + 'Rx' + str(receiver) + 'Refl' + str(reflectionOrder) + 'Trc' + traceIndex + '.csv'
    try:
        with open(filename) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
            xMpcCoordinate = []
            yMpcCoordinate = []
            zMpcCoordinate = []
            for row in readCSV:
                # Assign the coordinates
                xMpcCoordinate.append(row[::3])
                yMpcCoordinate.append(row[1::3])
                zMpcCoordinate.append(row[2::3])
        return xMpcCoordinate, yMpcCoordinate, zMpcCoordinate
    except FileNotFoundError:
        return None, None, None

#########################################
####    ENVIRONMENT COORDINATES      ####
#########################################
# Load the coordinates to display the environment topology (Provided by Q-D realization software)
def readEnvironmentCoordinates():
    fileName =  scenarioFolder + topologyFolder + topologyFilename
    x = []
    y = []
    z = []

    with open(fileName) as csvfile:
        topology = csv.reader(csvfile, delimiter=',',quoting=csv.QUOTE_NONNUMERIC)
        for row in topology:
            x.append(row[::3])
            y.append(row[1::3])
            z.append(row[2::3])
    return x,y,z


#########################################
#####    THROUGHPUT RESULTS      ########
#########################################
# Dictionnaries to handle the graphs for throughput
# Right now, it's only used for the  L-Room Scenario
THROUGHPUTMCS9_DIC = {} 
THROUGHPUTMCS12_DIC = {} 
THROUGHPUTMCS9_KEYSLIST = []
THROUGHPUTMCS9_VALUESLIST = []
THROUGHPUTMCS12_KEYSLIST = []
THROUGHPUTMCS12_VALUESLIST = []
# Load throughput values either for MCS9 or MCS12 (Provided by ns-3)
def readThroughputResults(mcs9):
    global THROUGHPUTMCS9_KEYSLIST
    global THROUGHPUTMCS9_VALUESLIST
    global THROUGHPUTMCS12_KEYSLIST
    global THROUGHPUTMCS12_VALUESLIST
    if mcs9 == True:
        filename = scenarioFolder + resultsFolder +  throughputMCS9File
    else:
        filename = scenarioFolder + resultsFolder +  throughputMCS12File
    try: 
        with open(filename) as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                if mcs9 == True:
                    THROUGHPUTMCS9_DIC[int(row['TRACE'])] = float(row['THROUGHPUT'])
                   
                else:
                    THROUGHPUTMCS12_DIC[int(row['TRACE'])] = float(row['THROUGHPUT'])
    except FileNotFoundError:
        print("Warning: Throughput file: ", filename, " does not exist")

    if mcs9 == True:
        THROUGHPUTMCS9_KEYSLIST= list(THROUGHPUTMCS9_DIC.keys())
        THROUGHPUTMCS9_VALUESLIST = list(THROUGHPUTMCS9_DIC.values())
       
    else:
        THROUGHPUTMCS12_KEYSLIST= list(THROUGHPUTMCS12_DIC.keys())
        THROUGHPUTMCS12_VALUESLIST = list(THROUGHPUTMCS12_DIC.values())
 
#########################################
########    SNR RESULTS      ############
#########################################
SNR_DIC = {} # Dictionnary containing the SNR values 
# Right now, it's only used for the L-Room Scenario
# Read the SNR results (Provided by ns-3)
def readSnrResults():
    filename = scenarioFolder + resultsFolder +   snrFile
    try: 
        with open(filename) as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                SNR_DIC[float(row['TRACE'])] = float(row['SNR'])
    except FileNotFoundError:
        print("Warning: Throughput file: ", filename, " does not exist")

   
#########################################
######    SLS PHASE RESULTS     #########
#########################################
SLS_DIC = {}
# Read the TxSS phase results (Provided by ns-3)
def readAllSLSResults():
    filename = scenarioFolder + resultsFolder +slsResultsFile
    try: 
        with open(filename) as f:
    
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                SLS_DIC[int(row['SRC_ID']), int(row['DST_ID']), int(row['TRACE_IDX'])] = int(row['SECTOR_ID'])
    except FileNotFoundError:
        print("Warning: SLS file: ", filename, " does not exist")
   
#########################################
##### ANTENNA ELEMENTS POSITION #########
#########################################
# Read the antenna elements position (Provided by Codebook Generator)
def readAntennaPosition():
    fileName = scenarioFolder + CodebookFolder + antennaArrangementFile
    with open(fileName, newline='') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
        antennaPosition = []
        for row in readCSV:
            antennaPosition.append(row)
    return antennaPosition
   
#########################################
#######    Nodes positions        #######
#########################################
# Read the file containing the position of each node in the topoloy (Provided by Q-D realization software)
def readAllPositionFile():
    global apCoordinates
    global staCoordinates
    global numberOfAp
    nbTraces = len(os.listdir(scenarioFolder + positionFolder))
    nodeCoordinatesTrace = []
    for i in range(nbTraces):
        colorNodes = []
        fileName = scenarioFolder + positionFolder +  positionFile + str(i) + ".csv" 
        with open(fileName, newline='') as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
            nodePosition = []
            for row in readCSV:
                nodePosition.append(row)
                # We want to color every nodes based on its ID
                colorNodes.append(readCSV.line_num)

        # Format the data as Mayavi is expecting them
        nodePositionZip = [list(i) for i in zip(*nodePosition)]
        nodePositionZip.append(colorNodes)
        nodeCoordinatesTrace.append(nodePositionZip)

    # Use numpy array to ease the split of STAs and APs coordinates
    npCoordinates = np.asarray(nodeCoordinatesTrace)

    # Split the APs and STAs coordinates
    apCoordinates = npCoordinates[::,::,:numberOfAp]
    staCoordinates = npCoordinates[::,::,numberOfAp:numberOfAp+numberOfSta]

    return apCoordinates,staCoordinates, nbTraces

#########################################
######    Antenna Patterns        #######
#########################################
# TODO Right now, we load the sme codebook for STAs and APs
# Load the codebook file (Provided by the Codebook Generator) to get the sectors directivity pattern 
def loadCodebook():
    fileName = scenarioFolder + CodebookFolder + CodebookFileName
    f = open(fileName)
    # The first line determines the number of phased antenna arrays within the device
    nbPhaseAntennaArrayNumber = int(f.readline())
    for antennaIndex in range (nbPhaseAntennaArrayNumber):
        
        # Read phased antenna array ID 
        antennaID = int(f.readline())
        # Read phased antenna array azimuth orientation degree
        azimuthOrientationDegree = int(f.readline())

        # Read phased antenna array elevation orientation degree
        elevationOrientationDegree = int(f.readline())

        # Read the number of antenna elements 
        elements =  int(f.readline())

        # Read the number of quantization bits for phase 
        phaseQuantizationBits = int(f.readline())
     
        # Read the number of quantization bits for amplitude */
        amplitudeQuantizationBits =  int(f.readline())

        # Read the directivity of a single antenna element
        singleElementDirectivity = np.empty((azimuthCardinality,elevationCardinality))
        
        for m in range (azimuthCardinality):
            line = f.readline().split(",")
            for n in range (elevationCardinality):
                singleElementDirectivity[m][n] = float(line[n]) 
       
        # Read the 3D steering vector of the antenna array
        steeringVector = np.zeros((azimuthCardinality,elevationCardinality,elements),dtype=np.complex128)
        for l in range (elements):
            for m in range (azimuthCardinality):
                line = f.readline().split(",")
                results = [float(i) for i in line]
                amp = results[::2]
                phaseDelay = results[1::2]
                
                for n in range(elevationCardinality):
                    steeringVector[m][n][l] = amp[n]*np.cos(phaseDelay[n]) + 1j*amp[n] * np.sin(phaseDelay[n])
         
        # Read Quasi-omni antenna weights 
        quasiOmniWeights = np.zeros(elements,dtype=np.complex128)
        line = f.readline().split(",")
        results = [float(i) for i in line]
        amp = results[::2]
        phaseDelay = results[1::2]
        for l in range (elements):
            polarSteering = cmath.polar(complex(amp[l],phaseDelay[l]))
            # It seems that there is an inconsitency between numpy and cmaths complex function
            quasiOmniWeights[l] = polarSteering[0] + 1j* polarSteering[1]

        # Read the number of sectors within this antenna array
        nbSectors = int(f.readline())
    
        for sector in range (nbSectors):
            sectorId =  int(f.readline())
            sectorType =  int(f.readline())
            sectorUsage =  int(f.readline())
            elementsWeights =  np.zeros(elements,dtype=np.complex128)
            line = f.readline().split(",")
            results = [float(i) for i in line]
            amp = results[::2]
            phaseDelay = results[1::2]
            for l in range (elements):
                # polarSteering = cmath.polar(complex(amp[l],phaseDelay[l]))
                # It seems that there is an inconsitency between numpy and cmaths complex function
                elementsWeights[l] = amp[l]*np.cos(phaseDelay[l]) + 1j*amp[l] * np.sin(phaseDelay[l])
            
            PATTERN_AP.append(calculateDirectivity(elementsWeights,steeringVector,singleElementDirectivity))
            PATTERN_STA.append(calculateDirectivity(elementsWeights,steeringVector,singleElementDirectivity))

# Calculate the directivity of each sector
def calculateDirectivity(weights,steeringVector,singleElementDirectivity):
    radiusFactor = 0.02 # Decide the effective size of the pattern when visualizing
    directivity = np.zeros((azimuthCardinality,elevationCardinality))
    # Try to speed-up the directivity computation by using just a slice of the azimuth to compute it
    directivityQuality = directivity[::quality,::1]
    steeringVectorQuality = steeringVector[::quality]
    singleElementDirectivityQuality = singleElementDirectivity[::quality] 
    azimuthQuality = directivityQuality.shape[0]
    elevationQuality = directivityQuality.shape[1]
    for m in range (azimuthQuality):
        for n in range (elevationQuality):
            value = complex(0)
            for j in range (len(weights)):
                value = value + weights[j] * steeringVectorQuality[m][n][j]
            value = value * singleElementDirectivityQuality[m][n]
            directivityQuality[m][n] = 10 * np.log10(abs(value)*abs(value))

    hi =  45
    lo =  -40
   
    # Clip the gain values
    powerDb = np.clip(directivityQuality,lo,hi)
    #Adjust for negative values
    if (lo < 0):
        for i in range(len(powerDb)):
            powerDb[i] -= lo

    azimuthAnglesWrapped = np.linspace(math.radians(0),math.radians(360),num=azimuthQuality)
    elevationAnglesWrapped =  np.linspace(math.radians(0),math.radians(180),num=elevationQuality)

    X, Y = np.meshgrid(elevationAnglesWrapped, azimuthAnglesWrapped)
    x = radiusFactor*powerDb*sin(X)*cos(Y)
    y = radiusFactor*powerDb*sin(X)*sin(Y)
    z = radiusFactor*powerDb*cos(X)

    return x,y,z,powerDb

######################################################################################################
#############                     END FUNCTIONS TO LOAD INPUT FILES                     ##############
######################################################################################################


######################################################################################################
###########                   FUNCTION TO IMPLEMENT THE VISUALIZATION                     ############
######################################################################################################

# Function to make an object visible
def makeVisible(anObject):
    anObject.actor.actor.visibility = True

# Function to make an object invisible
def makeInvisible(anObject):
    anObject.actor.actor.visibility = False

# Force the render (inspired by https://python.hotexamples.com/examples/pyface.api/GUI/set_busy/python-gui-set_busy-method-examples.html)
def force_render():
    _gui = GUI()
    orig_val = _gui.busy
    _gui.set_busy(busy=True)
    _gui.set_busy(busy=orig_val)
    _gui.process_events()

# Function to build a triangular mesh (used to build the environment)
def buildTriangle(triangleTexture, trianglesLimit,x,y,z, currentFigure, frontFaceCulling, nameComponent):
    textureFile = tvtk.JPEGReader()
    textureFile.file_name=triangleTexture + ".jpg" #any jpeg file
    triangleTexture = tvtk.Texture(input_connection=textureFile.output_port, interpolate=0)
    triangles = [(i*3, (i*3)+1, (i*3)+2 ) for i in range(trianglesLimit[0], trianglesLimit[1])]
    environmentTriangle = mlab.triangular_mesh(x, y, z,triangles,  colormap=colorMapCustom, color = (1,1,1),  representation = 'surface', opacity=1, figure = currentFigure, name = nameComponent)
    environmentTriangle.actor.enable_texture = True
    environmentTriangle.actor.tcoord_generator_mode = 'plane'
    environmentTriangle.actor.actor.texture = triangleTexture   
    environmentTriangle.actor.property.frontface_culling = frontFaceCulling   
    return environmentTriangle
    
# Function to display the antenna pattern directivity of the sector the AP uses to communicate with its associated STA 'node' 
def displayCurrentPatternAp(traceIndex,node):
    txssSector = -1
    realApIndex = 0

    # Find the sector used by the AP to transmit to the STA
    for apId in range(1,numberOfAp+1):
        if traceIndex == 0:
            if (apId,node, traceIndex) in SLS_DIC: 
                txssSector = SLS_DIC[(apId,node, traceIndex)]
                realApIndex = apId-1
        else:
            for i in range (traceIndex, -1 , -1):
                if (apId,node, i) in SLS_DIC: 
                    txssSector = SLS_DIC[(apId,node, i)]
                    realApIndex = apId -1
                    break

    if txssSector != -1:
        txssSector = txssSector-1 
        # The TxSS sector is known
        # Check if mobility was used in the scenario to correctly fill the position Index
        indexPosition = 0
        if nbTraces > 1:
            indexPosition = traceIndex
        XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodes = apCoordinates[indexPosition]
        
        # Update the antenna pattern position to the AP position
        ANTENNA_PATTERN_AP_VIEW2[txssSector].actor.actor.position = (XAPCOORDINATES[realApIndex],YAPCOORDINATES[realApIndex],ZAPCOORDINATES[realApIndex])

        # Color the antenna pattern
        ANTENNA_PATTERN_AP_VIEW2[txssSector].module_manager.scalar_lut_manager.data_range = [
                np.amin(PATTERN_AP[txssSector][3]), np.amax(PATTERN_AP[txssSector][3])] 
                
        # We want to hide the other AP antenna patterns
        for patternId in range(len(ANTENNA_PATTERN_AP_VIEW2)):
            makeInvisible(ANTENNA_PATTERN_AP_VIEW2[patternId])

        # and display the one corresponding to the AP TxSS
        makeVisible(ANTENNA_PATTERN_AP_VIEW2[txssSector]) 
        return realApIndex

# Function to display the antenna pattern directivity of the sector the STA 'node' uses to communicate with the AP it's associated with
def displayCurrentPatternSta(traceIndex,node):
    txssSector = -1
    for apId in range(1,numberOfAp+1):
        if traceIndex == 0:
            if (node,apId, traceIndex) in SLS_DIC:
                txssSector =  txssSector = SLS_DIC[(node,apId, traceIndex)]
        else:
            for i in range (traceIndex, -1 , -1):
                if (node,apId, i) in SLS_DIC:
                    txssSector = SLS_DIC[(node,apId, i)]
                    break
    
    if txssSector != -1:
        indexPosition = 0
        txssSector = txssSector-1
        if nbTraces > 1:
            indexPosition = traceIndex
        XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodes = staCoordinates[indexPosition]   
        realStaIndex = node - 1 - numberOfAp
        # Update the antenna pattern position to the STA current position
        ANTENNA_PATTERN_STA_VIEW2[txssSector].actor.actor.position = (XSTACOORDINATES[realStaIndex],YSTACOORDINATES[realStaIndex],ZSTACOORDINATES[realStaIndex])

        # Color the antenna pattern
        ANTENNA_PATTERN_STA_VIEW2[txssSector].module_manager.scalar_lut_manager.data_range = [
                np.amin(PATTERN_STA[txssSector][3]), np.amax(PATTERN_AP[txssSector][3])]       
    
        # We want to hide the other STA antenna patterns
        for patternId in range(len(ANTENNA_PATTERN_STA_VIEW2)):
            makeInvisible(ANTENNA_PATTERN_STA_VIEW2[patternId])
              
        # and display the one corresponding to the STA TxSS
        makeVisible(ANTENNA_PATTERN_STA_VIEW2[txssSector])
        return txssSector

# Display the AP antenna Geometry
def displayCurrentGeometryAP(traceIndex,apId):
    # Check if mobility was used in the scenario to correctly fill the position Index
    indexPosition = 0
    if nbTraces > 1:
        indexPosition = traceIndex

    # Update the coordinates of each element with the coordinates of the AP to center the antenna element on the AP
    XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodes = apCoordinates[indexPosition]
    # # Update the Geometry elements
    ANTENNA_GEOMETRY_AP_VIEW1[apId].actor.actor.position = (XAPCOORDINATES[apId],YAPCOORDINATES[apId],ZAPCOORDINATES[apId])

# Display the STA antenna Geometry
def displayCurrentGeometrySTA(traceIndex, node):
    # Check if mobility was used in the scenario to correctly fill the position Index
    indexPosition = 0
    if nbTraces > 1:
        indexPosition = traceIndex

    # Update the coordinates of each element with the coordinates of the STA
    XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodes = staCoordinates[indexPosition]
    realStaIndex = node - 1 - numberOfAp
    ANTENNA_GEOMETRY_STA_VIEW1[realStaIndex].actor.actor.position = (XSTACOORDINATES[realStaIndex],YSTACOORDINATES[realStaIndex],ZSTACOORDINATES[realStaIndex])
              

class Visualization(HasTraits):
    loadInput() # Parse input files and command-line paramter
   
    engine1 = Instance(Engine, args=())
   
    view1 = Instance(MlabSceneModel, ())
    view2 = Instance(MlabSceneModel, ())
    
    traceIndex = Range(0, nbTraces-1, 0) # File the trace index slider with the right number of traces

    # Initialize the STA enum based on the number of STAs
    STA = []
    STA.append("None")
    for i in range(numberOfAp,numberOfSta+numberOfAp):
        STA.append(i+1)

    STA_Training = Enum(*STA)

  
    ################################################################
    #################   INITIALIZE THE VIEWS      ##################
    ################################################################
    # By default, the environment geometry is displayed on view1 and view2
    # View1 displaysall the nodes (STAs+APs) and their antenna elements position
    # We are however computing everything a-priori all the antenna pattern for View 2
    # It takes more time to launch but is faster to update
    def __init__(self):
        # Do not forget to call the parent's __init__
        global AP_VIEW1, STA_VIEW1
        sceneLab = MlabSceneModel(engine=self.engine1) 

        # Set the background color for view1 and view2
        self.view2.background = (0.2901960784313726, 0.2901960784313726, 0.2901960784313726)
        self.view1.background = (0.2901960784313726, 0.2901960784313726, 0.2901960784313726)

        ################################################################
        #################            VIEW 1           ##################
        ################################################################

        ########################################################################################
        #############################    Display the topology    ###############################
        ########################################################################################

        # Create STA + AP (they are represented by 3D points)

        # Get the AP nodes coordinates
        XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodes = apCoordinates[self.traceIndex]
        # In order to properly manage the label, convert the nodes color to int
        colorNodesInt = colorNodes.astype(int)
        # Draw the APs nodes
        AP_VIEW1 = mlab.points3d(XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES,
                                                colorNodesInt, scale_factor=0.4, scale_mode="none", figure=self.view1.mayavi_scene,name = "AP_Nodes", reset_zoom = False)

        # Add a label for every AP
        labels = Labels()
        vtk_data_source = AP_VIEW1
        self.engine1.add_filter(labels, vtk_data_source)
        labels.mapper.label_format = ("AP %d")
        labels.mapper.label_mode = ('label_field_data')

        # Get the STA nodes coordinates
        XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodes = staCoordinates[self.traceIndex]
        # In order to properly manage the label, convert the nodes color to int
        colorNodesInt = colorNodes.astype(int)
        # Draw the STAs nodes
        STA_VIEW1 = mlab.points3d(XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES,
                                                colorNodesInt, scale_factor=0.4, scale_mode="none", figure=self.view1.mayavi_scene, name = "STA_Nodes",reset_zoom = False)

        # Add a label for every STA
        labels = Labels()
        vtk_data_source = STA_VIEW1
        self.engine1.add_filter(labels, vtk_data_source)
        labels.mapper.label_format = ("STA %d")
        labels.mapper.label_mode = ('label_field_data')
       
        ########################################################################################
        ###########################    Display antenna elements  ###############################
        ########################################################################################
        for i in range(numberOfAp):
            ANTENNA_GEOMETRY_AP_VIEW1.append(mlab.points3d(apAntennaGeometry[0][0],apAntennaGeometry[0][1], apAntennaGeometry[0][2],
                                                 scale_factor=0.1, scale_mode="none", figure=self.view1.mayavi_scene, name = "AP AntennaPosition",reset_zoom = False))

        # Display AP antenna elements arrangement
        for i in range(numberOfAp):
            displayCurrentGeometryAP(self.traceIndex,i)
        for i in range(numberOfSta):
            ANTENNA_GEOMETRY_STA_VIEW1.append(mlab.points3d(staAntennaGeometry[0][0],staAntennaGeometry[0][1], staAntennaGeometry[0][2],
                                                 scale_factor=0.1, scale_mode="none", figure=self.view1.mayavi_scene, name = "STA AntennaPosition",reset_zoom = False))
       
        # Display STA antenna elements arrangement
        for i in range(numberOfSta):
            displayCurrentGeometrySTA(self.traceIndex,i+2)

            
        ########################################################################################
        #######################    Display environment geometry  ###############################
        ########################################################################################

        
        x, y, z =  readEnvironmentCoordinates()
        for i in range(0,len(x),2):
            buildTriangle("hardfloor", [i,i+2],x, y, z,self.view1.mayavi_scene, True, "XWALL1") # Build the environment  for view1
            buildTriangle("hardfloor", [i,i+2],x, y, z,self.view2.mayavi_scene, True, "XWALL1")# Build the environment  for view2
    
        ########################################################################################
        ##########################    AP All TX Antenna Pattern    #############################
        ########################################################################################    

        # Create all the APs Tx Antenna Pattern
        colorValue = [0, 1]
        for antennaId in range(nbAntennaAp):
            for sectorId in range(nbSectorsPerAntenna):
                xAntennaPattern = PATTERN_AP[sectorId+antennaId*nbSectorsPerAntenna][0]
                yAntennaPattern = PATTERN_AP[sectorId+antennaId*nbSectorsPerAntenna][1]
                zAntennaPattern = PATTERN_AP[sectorId+antennaId*nbSectorsPerAntenna][2]
                colorAntennaPattern = PATTERN_AP[sectorId+antennaId*nbSectorsPerAntenna][3]
                ANTENNA_PATTERN_AP_VIEW2.append(mlab.mesh(xAntennaPattern, yAntennaPattern, zAntennaPattern,  vmin=min(colorValue), vmax=max(colorValue), tube_radius=0.025, figure=self.view2.mayavi_scene, name = "Antenna Pattern AP - Sector: "+str(sectorId+1) + " Antenna: " + str(antennaId+1), reset_zoom=False, scalars=colorAntennaPattern)
                                                )
                # Hide the antenna pattern
                ANTENNA_PATTERN_AP_VIEW2[sectorId+antennaId*nbSectorsPerAntenna].actor.actor.visibility = False
        ########################################################################################
        ##########################    AP STA TX Antenna Pattern    #############################
        ########################################################################################

        # Create all the STAs Tx Antenna Pattern
        for antennaId in range(nbAntennaSta):
            for sectorId in range(nbSectorsPerAntenna):
                xAntennaPattern = PATTERN_STA[sectorId+antennaId*nbSectorsPerAntenna][0]
                yAntennaPattern = PATTERN_STA[sectorId+antennaId*nbSectorsPerAntenna][1]
                zAntennaPattern = PATTERN_STA[sectorId+antennaId*nbSectorsPerAntenna][2]
                colorAntennaPattern = PATTERN_AP[sectorId+antennaId*nbSectorsPerAntenna][3] 
                ANTENNA_PATTERN_STA_VIEW2.append(mlab.mesh(xAntennaPattern, yAntennaPattern, zAntennaPattern,  vmin=min(colorValue), vmax=max(colorValue), tube_radius=0.025, figure=self.view2.mayavi_scene, name = "Antenna Pattern STA - Sector: "+str(sectorId+1) + " Antenna: " + str(antennaId+1), reset_zoom = False, scalars=colorAntennaPattern)
                                                )
                # Hide the antenna pattern
                ANTENNA_PATTERN_STA_VIEW2[sectorId+antennaId*nbSectorsPerAntenna].actor.actor.visibility = False

        HasTraits.__init__(self)


           
    ################################################################
    #            UPDATE WHEN THE SELECTED STA CHANGES              #
    ################################################################
    # VIEW2 is the only one to need to be updated as VIEW 1 displays all the nodes (STAs+APs)
    @on_trait_change('STA_Training')
    def update_SECOND_VIEW(self):
        global ANTENNA_PATTERN_AP_VIEW2
        global ANTENNA_PATTERN_STA_VIEW2
        global MPC
              
        if self.STA_Training == "None":
            # Remove everything from VIEW2
            
            # Hide the MPC
            for mpcIndex in range(len(MPC)):
                makeInvisible(MPC[mpcIndex])

            # Hide the AP Antenna Pattern
            for i in range(len(ANTENNA_PATTERN_AP_VIEW2)):
                makeInvisible(ANTENNA_PATTERN_AP_VIEW2[i])
            # Hide the STA Antenna Pattern
            for i in range(len(ANTENNA_PATTERN_STA_VIEW2)):
                makeInvisible(ANTENNA_PATTERN_STA_VIEW2[i])
        else:
            # A new STA is selected, update the display
            #####################################################
            ###########   Display STA and AP Txss   #############
            #####################################################
            indexPosition = 0
            if nbTraces > 1:
                indexPosition = self.traceIndex
            XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodesAp = apCoordinates[indexPosition]
            XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodesSta = staCoordinates[indexPosition]
            apIndex = displayCurrentPatternAp(self.traceIndex,self.STA_Training) # Display AP pattern
            displayCurrentPatternSta(self.traceIndex,self.STA_Training) # Display STA pattern

            #####################################################
            #############      Display MPCs         #############  
            #####################################################      
            for i in range (len(MPC)):
                MPC[i].remove()
            MPC[:] = [] 
            nbMpcToDraw = 0
            if len(MPC) == 0:

                # If the MPCs have never been displayed, create and display them
                for reflectionOrder in range(0, desiredReflectionOrder+1):
                    # Get the coordinates to draw the MPCs depending on the reflection order
                    xMpcCoordinate, yMpcCoordinate, zMpcCoordinate = readMPCCoordinateFile(
                        apIndex, self.STA_Training-1,reflectionOrder,str(indexPosition))
                    if reflectionOrder == 0:
                        # Display Direct Path
                        MPC.append(mlab.plot3d(xMpcCoordinate[0], yMpcCoordinate[0],
                                                zMpcCoordinate[0], tube_radius=0.1, colormap=colorMapCustom,color = (0,0,0), vmin=0, vmax=1,name = "Direct Path",  reset_zoom = False))
               
                        nbMpcToDraw = nbMpcToDraw +1
                    if reflectionOrder > 0:
                        # Draw the MPCs corresponding to the reflection order
                        for i in range(len(xMpcCoordinate)):
                            xs, ys, zs = xMpcCoordinate[i], yMpcCoordinate[i], zMpcCoordinate[i]
                            MPC.append(mlab.plot3d(
                                xs, ys, zs, tube_radius=0.025, colormap=colorMapCustom, color = (0,0,0), vmin=0, vmax=1, reset_zoom = False))
                            nbMpcToDraw = nbMpcToDraw +1
                        # Build 3 fake MPCS
                        MPC.append(mlab.plot3d(
                                xs, ys, zs, tube_radius=0.025, colormap=colorMapCustom, color = (0,0,0), vmin=0, vmax=1, reset_zoom = False))
                        MPC.append(mlab.plot3d(
                                xs, ys, zs, tube_radius=0.025, colormap=colorMapCustom, color = (0,0,0), vmin=0, vmax=1, reset_zoom = False))
                        MPC.append(mlab.plot3d(
                                xs, ys, zs, tube_radius=0.025, colormap=colorMapCustom, color = (0,0,0), vmin=0, vmax=1, reset_zoom = False))
   
            for mpcIndex in range(len(MPC),9):
                makeInvisible(MPC[mpcIndex]) 

        force_render()
        self.view1.render()
        self.view2.render()
        
    # Update the VIEW2 when the trace Index is changed
    @on_trait_change('traceIndex')
    def update_SECTOR(self):
        global ANTENNA_PATTERN_AP_VIEW2
        global ANTENNA_PATTERN_STA_VIEW2
        global MAPPING_RXSS_TO_SECTOR
        global MPC
        global THROUGHPUTMCS9_KEYSLIST , THROUGHPUTMCS9_VALUESLIST, THROUGHPUTMCS12_KEYSLIST, THROUGHPUTMCS12_VALUESLIST, graphData


        ################################################################
        #################       UPDATE VIEW 1         ##################
        ################################################################
        # Update the APs and STAs positions
        # Check if mobility was used in the scenario to correctly fill the position Index
        indexPosition = 0
        if nbTraces > 1:
            indexPosition = self.traceIndex
        # Update AP coordinates for the current trace
        XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodes = apCoordinates[indexPosition]
        XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodes = staCoordinates[indexPosition]
        # Update STA coordinates for the current trace
        STA_VIEW1.mlab_source.trait_set(x=XSTACOORDINATES, y = YSTACOORDINATES, z = ZSTACOORDINATES)
        AP_VIEW1.mlab_source.trait_set(x=XAPCOORDINATES, y = YAPCOORDINATES, z = ZAPCOORDINATES)



        # Update the antenna arrangement
        # Update APs antenna elements arrangement
        for i in range(numberOfAp):
            displayCurrentGeometryAP(self.traceIndex,i)
        # Update STAs antenna elements arrangement
        for i in range(numberOfSta):
            displayCurrentGeometrySTA(self.traceIndex,i+2)
      

        ################################################################
        #################       UPDATE VIEW 2         ##################
        ################################################################
      
        if self.STA_Training  != "None":
            # Check if mobility was used in the scenario to correctly fill the position Index
            indexPosition = 0
            if nbTraces > 1:
                indexPosition = self.traceIndex
           
            #####################################################
            ##########   Update STA and AP TxSS      ###########
            #####################################################
            # Grab the coordinates where to display AP and STA antenna pattern
            XAPCOORDINATES, YAPCOORDINATES, ZAPCOORDINATES, colorNodes = apCoordinates[indexPosition]
            XSTACOORDINATES, YSTACOORDINATES, ZSTACOORDINATES, colorNodes = staCoordinates[indexPosition]
            apIndex = displayCurrentPatternAp(self.traceIndex,self.STA_Training)  # Update AP Tx pattern
            displayCurrentPatternSta(self.traceIndex,self.STA_Training)  # Update STA Tx pattern
           
            
            #####################################################
            ################   Update MPCs      #################
            #####################################################
            hideDirectPath = False
            hideAllMpc = False
            if len(MPC) > 0:
                currentMPCIndex = 0
                # If the MPCs have never been displayed, create and display them
                for reflectionOrder in range(0, desiredReflectionOrder+1):
                    # Get the coordinates to draw the MPCs depending on the reflection order
                    xMpcCoordinate, yMpcCoordinate, zMpcCoordinate = readMPCCoordinateFile(
                        apIndex, self.STA_Training-1,reflectionOrder,str(indexPosition))
                    # Draw the MPC only if there is any
                    if xMpcCoordinate != None :
                        if reflectionOrder == 0:
                        # Direct Path
                            MPC[currentMPCIndex].mlab_source.trait_set(x=xMpcCoordinate[0], y=yMpcCoordinate[0], z=zMpcCoordinate[0])
                            currentMPCIndex = currentMPCIndex + 1
                        if reflectionOrder > 0:
                            # Draw the MPCs corresponding to the reflections
                            for i in range(len(xMpcCoordinate)):
                                xs, ys, zs = xMpcCoordinate[i], yMpcCoordinate[i], zMpcCoordinate[i]
                                MPC[currentMPCIndex].mlab_source.trait_set(x=xs, y=ys, z=zs)
                                currentMPCIndex = currentMPCIndex + 1
                    else:
                        currentMPCIndex = currentMPCIndex + 1
                        if reflectionOrder == 0:
                            hideDirectPath = True
                        else:
                            hideAllMpc = True
                            
            # Display the MPCs that must be displayed
            for mpcIndex in range(currentMPCIndex,len(MPC)):
                makeInvisible(MPC[mpcIndex])
            for mpcIndex in range(0,currentMPCIndex):
                if   hideDirectPath == True and mpcIndex == 0:
                    makeInvisible(MPC[mpcIndex])
                else:
                    makeVisible(MPC[mpcIndex])
            if hideAllMpc == True:
                for mpcIndex in range(0,currentMPCIndex):
                    makeInvisible(MPC[mpcIndex])

            
            #####################################################
            ###############   Update Graph      #################
            #####################################################
            if graphData:
            # Display MCS9 and MCS12 throughput
                xTime = []
                yMCS9 = []
                yMCS12 = []
                for i in range(0,self.traceIndex):
                    xTime.append(THROUGHPUTMCS9_KEYSLIST[i]/10)
                    yMCS9.append(THROUGHPUTMCS9_VALUESLIST[i])
                    yMCS12.append(THROUGHPUTMCS12_VALUESLIST[i])
                p2.setGeometry(p1.vb.sceneBoundingRect())
                p2.linkedViewChanged(p1.vb, p2.XAxis)
                p1.plot(xTime, yMCS9, pen=QPen(QColor(255, 0, 0)), clear = True)
                p1.plot(xTime, yMCS12,  pen=QPen(QColor(0, 153, 255)))

                # Add SNR values
                xSNR = []
                ySNR = []
                for key in sorted(SNR_DIC):
                    if key <= self.traceIndex/10:
                        xSNR.append(key)
                        ySNR.append(SNR_DIC[key])
            
                p2.clear()
                p2.addItem(pg.PlotDataItem(xSNR, ySNR, pen=QPen(QColor(0, 204, 0)), clear = True, name = "SNR"))

            force_render()
            self.view1.render()
            self.view2.render()
        
    # the layout of the dialog created
    view = View(             
        HSplit(
                Group(Item('view1', editor=SceneEditor(scene_class=MayaviScene),
                           height=0, width=0, show_label=False, resizable=False,),
                      ),
                Group(
                    Item('view2',
                         editor=SceneEditor(scene_class=MayaviScene), height=800,
                         width=600, show_label=False, resizable=False,),
                    VGroup(
                        Item(name='STA_Training',
                             style='simple'),
                        
                    ),
                    show_labels=False,
                ),
                ),
                '_', 'traceIndex',
                resizable=False,
                )



################################################################################
# The QWidget containing the visualization, this is pure PyQt4 code.
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)


        # There is a bug apparently in Mayavi that has not been fixed yet
        # (See https://github.com/enthought/mayavi/issues/3)
        # These 3 lines just catch the warning displayed due to this bug
        # and avoid to output them - Mainly cause it's speeding-up the simulation
        output=vtk.vtkFileOutputWindow()
        output.SetFileName("log.txt")
        vtk.vtkOutputWindow().SetInstance(output)
        self.visualization = Visualization()

        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        #QtCore.pyqtRemoveInputHook()
        #import pdb ; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

if __name__ == "__main__":
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    container.setWindowTitle("mmWave Python Visualizer")
    layout = QtGui.QGridLayout(container)

    mayavi_widget = MayaviQWidget(container)
    mayavi_widget.resize(1500,1000)
    
    layout.addWidget(mayavi_widget, 0, 0)

    # Create the widget to display the graph only if we have data to display
    if graphData:
        pg.setConfigOption('background', (68, 68, 68))
        plotWidget = pg.PlotWidget()
        plotWidget.addLegend(offset=(30,750))
        p1 = plotWidget.plotItem
        p1.setLabels(left='Throughput (Mbps)')
        p1.getAxis('left').setPen((255, 255, 255))
        p1.setYRange(0, 4000, padding=0)
        p2 = pg.ViewBox()
        p2.setYRange(0, 60, padding=0)
        p1.showAxis('right')
        p1.scene().addItem(p2)
        p1.getAxis('right').linkToView(p2)
        p2.setXLink(p1)
        p1.getAxis('right').setLabel('SNR (dB)', color='#00cc00')
        p1.getAxis('bottom').setPen((255, 255, 255))
        p2.setGeometry(p1.vb.sceneBoundingRect())
        p1.plot(name = "Throughput MCS9", pen=QPen(QColor(255, 0, 0)))
        p1.plot(name = "Throughput MCS12",pen=QPen(QColor(0, 153, 255)))
        p1.plot(name = "SNR",pen=QPen(QColor(0, 204, 0)))
        ## need to re-update linked axes since this was called
        ## incorrectly while views had different shapes.
        ## (probably this should be handled in ViewBox.resizeEvent)
        p2.linkedViewChanged(p1.vb, p2.XAxis)
        layout.addWidget(plotWidget, 0, 1)  

    
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()

    # Start the main event loop.
    app.exec_()