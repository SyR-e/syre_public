import pickle
import sys
import os
import os.path
import math

currentFolder = os.path.dirname(os.path.realpath(__file__))

filepath = r'%s/temp/temp.pkl'%(currentFolder)
with open(filepath,'rb') as exportdata:
    simdata = pickle.load(exportdata)

sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64")
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64/PythonFiles/DesktopPlugin")
import ScriptEnv

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("%s"%(simdata["filename"][:-4]))
oDesign = oProject.SetActiveDesign("Maxwell2DDesign1")
#oDesign = oProject.SetActiveDesign("SOLID_ACLOSS")
oEditor = oDesign.SetActiveEditor("3D Modeler")

#oDesktop.AddMessage ("%s"%(geo["filename"][:-4]), "Maxwell2DDesign1", 0, "%s"%(simdata["iAmp"]), "")

############## Current assignment excitation per phases

oModule = oDesign.GetModule("BoundarySetup")
oModule.DeleteAllExcitations()
oModule.AssignWindingGroup(
	[
		"NAME:WG_Ph1",
		"Type:="		, "Current",
		"IsSolid:="		, False,
		"Current:="		, "%s*cos(%s*2*pi/60*%s*time+(%s))A"%(simdata["iAmp"],simdata["EvalSpeed"],simdata["p"],simdata["phi_init"]), 
		"Resistance:="		, "0ohm",
		"Inductance:="		, "0nH",
		"Voltage:="		, "0mV",
		"ParallelBranchesNum:="	, "1"
	])
oModule.AssignWindingGroup(
	[
		"NAME:WG_Ph2",
		"Type:="		, "Current",
		"IsSolid:="		, False,
		"Current:="		, "%s*cos(%s*2*pi/60*%s*time+(240)*pi/180+(%s))A"%(simdata["iAmp"],simdata["EvalSpeed"],simdata["p"],simdata["phi_init"]), 
		"Resistance:="		, "0ohm",
		"Inductance:="		, "0nH",
		"Voltage:="		, "0mV",
		"ParallelBranchesNum:="	, "1"
	])
oModule.AssignWindingGroup(
	[
		"NAME:WG_Ph3",
		"Type:="		, "Current",
		"IsSolid:="		, False,
		"Current:="		, "%s*cos(%s*2*pi/60*%s*time+(120)*pi/180+(%s))A"%(simdata["iAmp"],simdata["EvalSpeed"],simdata["p"],simdata["phi_init"]), 
		"Resistance:="		, "0ohm",
		"Inductance:="		, "0nH",
		"Voltage:="		, "0mV",
		"ParallelBranchesNum:="	, "1"
	])




############ lists name slots 
ii=0
slotconds=""
lineslotairgaps=""
slotairgaps=""
statplate="1_1"

for ii in range(simdata["q"]*3):  

	slotindex=2+ii*4
	lineindex=slotindex+1
	slot="1_%d" %(slotindex)
	line="1_%d"%(lineindex)
	lineairgap="1_%d,1_%d"%(lineindex+1,lineindex+2)
	slotconds = slotconds + "," + slot + "," + slot + "_Detach1"   ##list slot/cond, useful for detachment
	lineslotairgaps=lineslotairgaps+","+lineairgap
	slotairgaps=slotairgaps+","+statplate+"_Detach"+"%d"%(ii+1)

slotconds=slotconds[1:] #removing initial comma of the string
slotconds_array=slotconds.split(',') #vector containig the elements of the string separated by the comma


for jj in range(len(simdata["win_avv"][0])): #n slot per layer

	for ii in range(len(simdata["win_avv"])): # n layer slot
	
		if simdata["win_avv"][ii][jj]>0:
			oModule.AssignCoil(
				[
					"NAME: Coil_%s"%(slotconds_array[len(simdata["win_avv"])*jj+ii]),
					"Objects:="		, [slotconds_array[len(simdata["win_avv"])*jj+ii]],
					"Conductor number:="	, "%s"%(simdata["N_cond"]),
					"PolarityType:="	, "Positive"
				])
		else:
			oModule.AssignCoil(
				[
					"NAME: Coil_%s"%(slotconds_array[len(simdata["win_avv"])*jj+ii]),
					"Objects:="		, [slotconds_array[len(simdata["win_avv"])*jj+ii]],
					"Conductor number:="	, "%s"%(simdata["N_cond"]),
					"PolarityType:="	, "Negative"
				])

		phasenum=abs(int(simdata["win_avv"][ii][jj]))
		oModule.AddWindingCoils("WG_Ph%d"%(phasenum), ["Coil_%s"%(slotconds_array[len(simdata["win_avv"])*jj+ii])])

####### Magnets Temperature ###########
if simdata["n_PM"]!=0:
	PMs_temp=""
	temp=0
	for jj in range (2):
		for ii in range(simdata["nlay"]+simdata["radial_ribs_split"]):
			rotorindex=(3+5*simdata["q"]*3)*(1-jj)+ii+temp*jj
			
		for kk in range(int(simdata["n_PM"]/2)):
			rotorindex=rotorindex+1
			PM="2_%d" %(rotorindex)
			PMs_temp=PMs_temp + "," + PM + "," + "%fcel"%(simdata["tempPP"])

		temp=rotorindex+1
	PMs_temp=PMs_temp[1:]
	PMs_temp_array=PMs_temp.split(',')


	oDesign.SetObjectTemperature(
		[
			"NAME:TemperatureSettings",
			"IncludeTemperatureDependence:=", True,
			"EnableFeedback:="	, False,
			"Temperatures:="	, PMs_temp_array
		])

####### Setup Core Loss Part ######
if simdata['corelossflag']==1:
	#rotorindex=3+5*simdata["q"]*3+2*(simdata["nlay"]+simdata["radial_ribs_split"])+simdata["n_PM"]
	rotorindex=3+5*simdata["q"]*3+2*(simdata["nlay"]+simdata["radial_ribs_split"])+simdata["n_PM"]
	#rotorindex= 41;
	rotorplate="2_%d" %(rotorindex)
	oModule.SetCoreLoss([statplate,rotorplate], False)
	#if shaft!="ShaftAir":
	#   oModule.SetCoreLoss([shaft], False)

####### Setup Moving part ##########

oModule = oDesign.GetModule("ModelSetup")
MotionSetupName = oModule.GetMotionSetupNames ()

if MotionSetupName == ["MotionSetup1"]:
	oModule = oDesign.GetModule("ReportSetup")
	oModule.DeleteAllReports()
	oModule = oDesign.GetModule("ModelSetup")
	oModule.DeleteMotionSetup(["MotionSetup1"])

oModule.AssignBand(
	[
		"NAME:Data",
		"Move Type:="		, "Rotate",
		"Coordinate System:="	, "Global",
		"Axis:="		, "Z",
		"Is Positive:="		, True,
		"InitPos:="		, "0deg",   ## Initial Rotor Position
		"HasRotateLimit:="	, False,
		"NonCylindrical:="	, False,
		"Consider Mechanical Transient:=", False,
		"Angular Velocity:="	, "%srpm"%(simdata["EvalSpeed"]),  ## Rotor Speed
		"Objects:="		, ["Rotating_Band_out"]
	])


####### Analisys Setup ##########

if simdata["save_fields"]==0:
	savefieldstype= "None"
else:
	savefieldstype= "Every N Steps"


oModule = oDesign.GetModule("AnalysisSetup")
oModule.ResetAllToTimeZero()

SetupName = oModule.GetSetups()
#oDesktop.AddMessage ("%s"%(simdata["filename"][:-4]), "Maxwell2DDesign1", 0, "%s"%(SetupName), "")
if SetupName == ["Setup1"]:
	oModule.DeleteSetups(["Setup1"])
	oDesktop.ClearMessages("", "",2)

oModule.InsertSetup("Transient", 
	[
		"NAME:Setup1",
		"Enabled:="		, True,
		[
			"NAME:MeshLink",
			"ImportMesh:="		, False
		],
		"NonlinearSolverResidual:=", "0.001",
		"TimeIntegrationMethod:=", "BackwardEuler",
		"SmoothBHCurve:="	, False,
		"StopTime:="		, "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"]),
		"TimeStep:="		,  "60/(%s*%s*360/%s*%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"],simdata["nsim"]),
		"OutputError:="		, False,
		"UseControlProgram:="	, False,
		"ControlProgramName:="	, " ",
		"ControlProgramArg:="	, " ",
		"CallCtrlProgAfterLastStep:=", False,
		"FastReachSteadyState:=", False,
		"AutoDetectSteadyState:=", False,
		"IsGeneralTransient:="	, True,
		"IsHalfPeriodicTransient:=", False,
		"SaveFieldsType:="	, "%s"%(savefieldstype),
		"N Steps:="		, "%s"%(simdata["save_fields"]),
		"Steps From:="		, "0s",
		"Steps To:="		, "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"]),
		"CacheSaveKind:="	, "Count",
		"NumberSolveSteps:="	, 1,
		"RangeStart:="		, "0s",
		"RangeEnd:="		, "0.1s",
		"UseAdaptiveTimeStep:="	, False,
		"InitialTimeStep:="	, "0.002s",
		"MinTimeStep:="		, "0.001s",
		"MaxTimeStep:="		, "0.003s",
		"TimeStepErrTolerance:=", 0.0001
	])

############# Solutions Plot #############

#Torque
oModule = oDesign.GetModule("ReportSetup")

oModule.CreateReport("Torque Plot", "Transient", "Rectangular Plot", "Setup1 : Transient", 
	[
		"Domain:="		, "Sweep"
	], 
	[
		"Time:="		, ["All"]
	], 
	[
		"X Component:="		, "Time",
		"Y Component:="		, ["Moving1.Torque"]
	])
oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Torque Plot:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "NewtonMeter"]]]])

#Induced Voltages
oModule.CreateReport("Induced Voltages", "Transient", "Rectangular Plot", "Setup1 : Transient", 
	[
		"Domain:="		, "Sweep"
	], 
	[
		"Time:="		, ["All"]
	], 
	[
		"X Component:="		, "Time",
		"Y Component:="		, ["InducedVoltage(WG_Ph1)","InducedVoltage(WG_Ph2)","InducedVoltage(WG_Ph3)"]
	])
oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Induced Voltages:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "V"]]]])

#Input Currents
oModule.CreateReport("Input Currents", "Transient", "Rectangular Plot", "Setup1 : Transient", 
	[
		"Domain:="		, "Sweep"
	], 
	[
		"Time:="		, ["All"]
	], 
	[
		"X Component:="		, "Time",
		"Y Component:="		, ["InputCurrent(WG_Ph1)","InputCurrent(WG_Ph2)","InputCurrent(WG_Ph3)"]
	])
oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Input Currents:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "A"]]]])

#Flux
oModule.CreateReport("Fluxes Linkages", "Transient", "Rectangular Plot", "Setup1 : Transient", 
	[
		"Domain:="		, "Sweep"
	], 
	[
		"Time:="		, ["All"]
	], 
	[
		"X Component:="		, "Time",
		"Y Component:="		, ["FluxLinkage(WG_Ph1)","FluxLinkage(WG_Ph2)","FluxLinkage(WG_Ph3)"]
	])
oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Fluxes Linkages:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "Wb"]]]])

#Mechanical Position Theta
oModule.CreateReport("Mechanical Position", "Transient", "Rectangular Plot", "Setup1 : Transient", 
	[
		"Domain:="		, "Sweep"
	], 
	[
		"Time:="		, ["All"]
	], 
	[
		"X Component:="		, "Time",
		"Y Component:="		, ["Moving1.Position"]
	])
oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Mechanical Position:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "deg"]]]])

####### Core Loss
if simdata['corelossflag']==1:
	
	oModule.CreateReport("Core Loss", "Transient", "Rectangular Plot", "Setup1 : Transient", 
		[
			"Domain:="		, "Sweep"
		], 
		[
			"Time:="		, ["All"]
		], 
		[
			"X Component:="		, "Time",
			"Y Component:="		, ["CoreLoss","CoreLoss(%s)"%(statplate),"CoreLoss(%s)"%(rotorplate)]
		])
	oModule.AddTraceCharacteristics("Core Loss", "avg", [], ["Specified", "60/(%s*%s)s"%(simdata['EvalSpeed'],simdata['p']), "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"])])
	oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Core Loss:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "W"]]]])

	oModule.CreateReport("Hysteresis Loss", "Transient", "Rectangular Plot", "Setup1 : Transient", 
		[
			"Domain:="		, "Sweep"
		], 
		[
			"Time:="		, ["All"]
		], 
		[
			"X Component:="		, "Time",
			"Y Component:="		, ["HysteresisLoss","HysteresisLoss(%s)"%(statplate),"HysteresisLoss(%s)"%(rotorplate)]
		])	
	oModule.AddTraceCharacteristics("Hysteresis Loss", "avg", [], ["Specified", "60/(%s*%s)s"%(simdata['EvalSpeed'],simdata['p']), "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"])])
	oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Hysteresis Loss:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "W"]]]])

	oModule.CreateReport("Eddy Currents Loss", "Transient", "Rectangular Plot", "Setup1 : Transient", 
		[
			"Domain:="		, "Sweep"
		], 
		[
			"Time:="		, ["All"]
		], 
		[
			"X Component:="		, "Time",
			"Y Component:="		, ["EddyCurrentLoss","EddyCurrentLoss(%s)"%(statplate),"EddyCurrentLoss(%s)"%(rotorplate)]
		])
	oModule.AddTraceCharacteristics("Eddy Currents Loss", "avg", [], ["Specified", "60/(%s*%s)s"%(simdata['EvalSpeed'],simdata['p']), "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"])])
	oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Eddy Currents Loss:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "W"]]]])

	oModule.CreateReport("Excess Loss", "Transient", "Rectangular Plot", "Setup1 : Transient", 
		[
			"Domain:="		, "Sweep"
		], 
		[
			"Time:="		, ["All"]
		], 
		[
			"X Component:="		, "Time",
			"Y Component:="		, ["ExcessLoss","ExcessLoss(%s)"%(statplate),"ExcessLoss(%s)"%(rotorplate)]
		])
	oModule.AddTraceCharacteristics("Excess Loss", "avg", [], ["Specified", "60/(%s*%s)s"%(simdata['EvalSpeed'],simdata['p']), "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"])])
	oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","Excess Loss:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "W"]]]])

oModule.CreateReport("SolidLoss", "Transient", "Rectangular Plot", "Setup1 : Transient", 
		[
			"Domain:="		, "Sweep"
		], 
		[
			"Time:="		, ["All"]
		], 
		[
			"X Component:="		, "Time",
			"Y Component:="		, ["SolidLoss"]
		])
	#oModule.AddTraceCharacteristics("SolidLoss", "avg", [], ["Specified", "60/(%s*%s)s"%(simdata['EvalSpeed'],simdata['p']), "60/(%s*%s*360/%s) s"%(simdata['EvalSpeed'],simdata['p'],simdata["xdeg"])])
	#oModule.ChangeProperty(["NAME:AllTabs",["NAME:Scaling",["NAME:PropServers","SolidLoss:AxisY1"],["NAME:ChangedProps",["NAME:Auto Units",	"Value:=", False],["NAME:Units","Value:=", "W"]]]])


########### Initial Mesh Settings ################

oModule = oDesign.GetModule("MeshSetup")
oModule.InitialMeshSettings(
	[
		"NAME:MeshSettings",
		[
			"NAME:GlobalSurfApproximation",
			"CurvedSurfaceApproxChoice:=", "UseSlider",
			"SliderMeshSettings:="	, 5
		],
		[
			"NAME:GlobalModelRes",
			"UseAutoLength:="	, True
		],
		"MeshMethod:="		, "AnsoftClassic"
	])
########### Launch Analisys ##############

oDesign.AnalyzeAll()

########### Export Plots #################
resultspath="%s%s_results/FEA results/Ansys_%s_%sA"%(simdata["filepath"],simdata["filename"][:-4],simdata["eval_type"],simdata['iAmp'])
if not os.path.exists(resultspath):
	os.makedirs(resultspath)

oModule = oDesign.GetModule("ReportSetup")
oModule.ExportToFile("Torque Plot", "%s/TorquePlotData.csv"%(resultspath), False)
oModule.ExportToFile("Induced Voltages", "%s/VoltagesPlotData.csv"%(resultspath), False)
oModule.ExportToFile("Input Currents", "%s/CurrentsPlotData.csv"%(resultspath), False)
oModule.ExportToFile("Fluxes Linkages", "%s/FluxesPlotData.csv"%(resultspath), False)
oModule.ExportToFile("Mechanical Position", "%s/PositionData.csv"%(resultspath), False)

if simdata['corelossflag']==1:
	#first graph value
	oModule.ExportToFile("Core Loss", "%s/CoreLossData.csv"%(resultspath), False)
	oModule.ExportToFile("SolidLoss", "%s/MagnetLossData.csv"%(resultspath), False)
	#avg value
	#oModule.ExportTableToFile("Core Loss", "%s/CoreLossAvg.txt"%(resultspath), "Legend")
	#oModule.ExportTableToFile("Hysteresis Loss", "%s/HysteresisLossAvg.txt"%(resultspath), "Legend")
	#oModule.ExportTableToFile("Eddy Currents Loss", "%s/EddyCurrentsLossAvg.txt"%(resultspath), "Legend")
	#oModule.ExportTableToFile("Excess Loss", "%s/ExcessLossAvg.txt"%(resultspath), "Legend")


sys.exit()
