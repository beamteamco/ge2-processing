# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.13-2 replay file
# Internal Version: 2013_07_18-05.24.06 126428
# Run by software-dev on Tue Jan 19 13:30:24 2016
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=672.394470214844, 
    height=212.765670776367)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
o1 = session.openOdb(name='/home/software-dev/MicroTest_1.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/software-dev/MicroTest_1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       2
#: Number of Node Sets:          7
#: Number of Steps:              3
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=2, frame=0 )
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=1, frame=2 )
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
odb = session.mdbData['Model-1']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
odb = session.odbs['/home/software-dev/MicroTest_1.odb']
session.viewports['Viewport: 1'].setValues(displayedObject=odb)
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=2, frame=605 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=2, frame=605 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=2, frame=604 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=2, frame=603 )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='RF', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=CONTOURS_ON_DEF)
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='LE', outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, 'Max. Principal'), )
#: Warning: The selected Primary Variable is not available in the current frame for any elements in the current display group.
