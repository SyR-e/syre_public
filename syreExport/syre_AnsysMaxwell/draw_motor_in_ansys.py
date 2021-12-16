import sys
import pickle
import os
currentFolder = os.path.dirname(os.path.realpath(__file__))
sys.path.append(r"C:\\Program Files\\AnsysEM\\AnsysEM20.1\\Win64")
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64/PythonFiles/DesktopPlugin")
import math

############### export data from mat file (through pickle) ################
filepath = r'%s/temp/temp.pkl'%(currentFolder)
with open(filepath,'rb') as exportdata:
    geo,material = pickle.load(exportdata)

import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.SetActiveProject("%s"%(geo["filename"][:-4]))
oDesign = oProject.SetActiveDesign("Maxwell2DDesign1")
oEditor = oDesign.SetActiveEditor("3D Modeler")



############### import dxf by file directory #############
oDesktop.AddMessage ("%s"%(geo["filename"][:-4]), "Maxwell2DDesign1", 0, "%s%s"%(geo["filepath"],geo["filename"][:-4]), "")
oEditor.ImportDXF(
	[
		"NAME:options",
		"FileName:="		, "%s%s.dxf" %(geo["filepath"],geo["filename"][:-4]),
		"Scale:="		, 0.001,
		"AutoDetectClosed:="	, True,
		"SelfStitch:="		, True,
		"DefeatureGeometry:="	, True,
		"DefeatureDistance:="	, 0.0001,
		"RoundCoordinates:="	, False,
		"RoundNumDigits:="	, 4,
		"SelfStitchTolerance:="	, 0.0001,
		"WritePolyWithWidthAsFilledPoly:=", True,
		"ImportMethod:="	, 1,
		"2DSheetBodies:="	, True,
		[
			"NAME:LayerInfo",
			[
				"NAME:0",
				"source:="		, "0",
				"display_source:="	, "0",
				"import:="		, True,
				"dest:="		, "0",
				"dest_selected:="	, False,
				"layer_type:="		, "signal"
			],
			[
				"NAME:1",
				"source:="		, "1",
				"display_source:="	, "1",
				"import:="		, True,
				"dest:="		, "1",
				"dest_selected:="	, False,
				"layer_type:="		, "signal"
			],
			[
				"NAME:2",
				"source:="		, "2",
				"display_source:="	, "2",
				"import:="		, True,
				"dest:="		, "2",
				"dest_selected:="	, False,
				"layer_type:="		, "signal"
			]
		]
	])

oEditor.FitAll()
####################### detach slot  ################
ii=0
slotconds=""
lineslotairgaps=""
slotairgaps=""
statplate="1_1"

for ii in range(geo["q"]*3):    
	slotindex=2+ii*4
	lineindex=slotindex+1
	slot="1_%d" %(slotindex)
	line="1_%d"%(lineindex)
	lineairgap="1_%d,1_%d"%(lineindex+1,lineindex+2)
	slotconds = slotconds + "," + slot + "," + slot + "_Detach1"   ##list slot/cond, useful for detachment
	lineslotairgaps=lineslotairgaps+","+lineairgap
	slotairgaps=slotairgaps+","+statplate+"_Detach"+"%d"%(ii+1)


	########## Slot
	oEditor.Imprint(
		[
			"NAME:Selections",
			"Blank Parts:="		, slot,
			"Tool Parts:="		, line
		], 
		[
			"NAME:ImprintParameters",
			"KeepOriginals:="	, False
		])

	facesIDs=oEditor.GetFaceIDs(slot)
	
	faceID=int(facesIDs[1])
	oEditor.DetachFaces(
		[
			"NAME:Selections",
			"Selections:="		, slot,
			"NewPartsModelFlag:="	, "Model"
		], 
		[
			"NAME:Parameters",
			[
				"NAME:DetachFacesToParameters",
				"FacesToDetach:="	, [faceID]
			]
		])

lineslotairgaps=lineslotairgaps[1:]

slotconds=slotconds[1:] #removing initial comma of the string
slotconds_array=slotconds.split(',') #vector containig the elements of the string separated by the comma

slotairgaps=slotairgaps[1:]
slotairgaps_array=slotairgaps.split(',')

############ detach slot from stator's iron ###############
oEditor.Subtract(
	[
		"NAME:Selections",
		"Blank Parts:="		, statplate,
		"Tool Parts:="		, slotconds
	], 
	[
		"NAME:SubtractParameters",
		"KeepOriginals:="	, True
	])
########### detach slot airgap from stator's iron


oEditor.Imprint(
	[
		"NAME:Selections",
		"Blank Parts:="		, statplate,
		"Tool Parts:="		, lineslotairgaps
	], 
	[
		"NAME:ImprintParameters",
		"KeepOriginals:="	, False
	])

facesIDs=oEditor.GetFaceIDs(statplate)
facesIDs=facesIDs[:-1] #remove last term
facesIDs=map(int,facesIDs) #array of string -> integer array

oEditor.DetachFaces(
	[
		"NAME:Selections",
		"Selections:="		, statplate,
		"NewPartsModelFlag:="	, "Model"
	], 
	[
		"NAME:Parameters",
		[
			"NAME:DetachFacesToParameters",
			"FacesToDetach:="	, facesIDs
		]
	])



	############ Subract barrier from rotor's iron ###############
barriers=""
PMs=""
temp=0
for jj in range (2):
	for ii in range(geo["nlay"]+geo["radial_ribs_split"]):
		rotorindex=(3+5*geo["q"]*3)*(1-jj)+ii+temp*jj
		barrier="2_%d" %(rotorindex)
		barriers=barriers+","+barrier
		
		

	for kk in range(int(geo["n_PM"]/2)):
		rotorindex=rotorindex+1
		PM="2_%d" %(rotorindex)
		PMs=PMs+","+PM

	temp=rotorindex+1
#oDesktop.AddMessage ("%s"%(geo["filename"][:-4]), "Maxwell2DDesign1", 0, "%s"%(barriers), "")
barriers=barriers[1:] 
barriers_array = barriers.split(",")
PMs=PMs[1:]
PMs_array=PMs.split(',')

rotorplate="2_%d" %(temp)
shaft="2_%d" %(temp+1)

oEditor.Subtract(
	[
		"NAME:Selections",
		"Blank Parts:="		, rotorplate,
		"Tool Parts:="		, barriers
	], 
	[
		"NAME:SubtractParameters",
		"KeepOriginals:="	, False
	])

if geo["n_PM"]!=0:
	oEditor.Subtract(
		[
			"NAME:Selections",
			"Blank Parts:="		, rotorplate,
			"Tool Parts:="		, PMs
		], 
		[
			"NAME:SubtractParameters",
			"KeepOriginals:="	, True
		])

################### Material Assignment ####################
###Rotor
oEditor.ChangeProperty(
	[
		"NAME:AllTabs",
		[
			"NAME:Geometry3DAttributeTab",
			[
				"NAME:PropServers",rotorplate
			],
			[
				"NAME:ChangedProps",
				[
					"NAME:Solve Inside","Value:="		, True
				],
				[
					"NAME:Color","R:=", 140,"G:=", 160,	"B:=", 175
				],
				[
					"NAME:Material",
					"Value:="		, "\"%s\"" %(material["rotor"])
				]
			]
		]
	])



###Stator
oEditor.ChangeProperty(
	[
		"NAME:AllTabs",
		[
			"NAME:Geometry3DAttributeTab",
			[
				"NAME:PropServers",statplate
			],
			[
				"NAME:ChangedProps",
				[
					"NAME:Solve Inside","Value:="		, True
				],
				[
					"NAME:Color","R:=", 140,"G:=", 160,	"B:=", 175
				],
				[
					"NAME:Material",
					"Value:="		, "\"%s\"" %(material["stator"])
				]
			]
		]
	])
	
###Shaft
oEditor.ChangeProperty(
	[
		"NAME:AllTabs",
		[
			"NAME:Geometry3DAttributeTab",
			[
				"NAME:PropServers",shaft
			],
			[
				"NAME:ChangedProps",
				[
					"NAME:Solve Inside","Value:="		, True
				],
				[
					"NAME:Color","R:=", 140,"G:=", 160,	"B:=", 175
				],
				[
					"NAME:Material",
					"Value:="		, "\"%s\"" %(material["shaft"])
				]
			]
		]
	])

oEditor.SetWCS(["NAME:SetWCS Parameter","Working Coordinate System:=", "Global","RegionDepCSOk:=", False])
###Magnet
if geo["n_PM"]!=0:
	
	for ii in range(len(PMs_array)):
		oEditor.CreateRelativeCS(
		[
			"NAME:RelativeCSParameters",
			"Mode:="		, "Axis/Position",
			"OriginX:="		, "%smm"%(geo["PM_CS"][0][ii]),
			"OriginY:="		, "%smm"%(geo["PM_CS"][1][ii]),
			"OriginZ:="		, "0mm",
			"XAxisXvec:="		, "%smm"%(geo["PM_CS"][2][ii]),
			"XAxisYvec:="		, "%smm"%(geo["PM_CS"][3][ii]),
			"XAxisZvec:="		, "0mm",
			"YAxisXvec:="		, "%smm"%(-geo["PM_CS"][3][ii]),
			"YAxisYvec:="		, "%smm"%(geo["PM_CS"][2][ii]),
			"YAxisZvec:="		, "0mm"
		], 
		[
			"NAME:Attributes",
			"Name:="		, "RelativeCS%s"%(PMs_array[ii])
		])

		oEditor.ChangeProperty(
			[
				"NAME:AllTabs",
				[
					"NAME:Geometry3DAttributeTab",
					[
						"NAME:PropServers",PMs_array[ii]
					],
					[
						"NAME:ChangedProps",
						[
							"NAME:Solve Inside","Value:="		, True
						],
						[
							"NAME:Color","R:=", 50,"G:=", 50,	"B:=", 50
						],
						[
							"NAME:Material",
							"Value:="		, "\"%s\"" %(material["magnet"])
						],
						[
							"NAME:Orientation",
							"Value:="		, "RelativeCS%s"%(PMs_array[ii])
						]
					]
				]
			])
		oEditor.SetWCS(["NAME:SetWCS Parameter","Working Coordinate System:=", "Global","RegionDepCSOk:=", False])
###Conductor
for ii in range(len(slotconds_array)): ##3 fasi
	oEditor.ChangeProperty(
		[
			"NAME:AllTabs",
			[
				"NAME:Geometry3DAttributeTab",
				[
					"NAME:PropServers",slotconds_array[ii]
				],
				[
					"NAME:ChangedProps",
					[
						"NAME:Solve Inside","Value:="		, True
					],
					[
						"NAME:Color","R:=", 220,"G:=", 165,	"B:=", 30
					],
					[
						"NAME:Material",
						"Value:="		, "\"%s\"" %(material["slotcond"])
					]
				]
			]
	])

# barrier removed because air background, programm more stable
# ###Aria: Barriere di flusso
# for ii in range(len(barriers_array)):
# 	oEditor.ChangeProperty(
# 		[
# 			"NAME:AllTabs",
# 			[
# 				"NAME:Geometry3DAttributeTab",
# 				[
# 					"NAME:PropServers",barriers_array[ii] 
# 				],
# 				[
# 					"NAME:ChangedProps",
# 					[
# 						"NAME:Solve Inside","Value:="		, True
# 					],
# 					[
# 						"NAME:Color","R:=", 80,"G:=", 150,	"B:=", 250
# 					],
# 					[
# 						"NAME:Material",
# 						"Value:="		, "\"Air\""
# 					]
# 				]
# 			]
# 		])

###Aria: Traferro di cava
for ii in range(len(slotairgaps_array)):
	oEditor.ChangeProperty(
		[
			"NAME:AllTabs",
			[
				"NAME:Geometry3DAttributeTab",
				[
					"NAME:PropServers",slotairgaps_array[ii]  
				],
				[
					"NAME:ChangedProps",
					[
						"NAME:Solve Inside","Value:="		, True
					],
					[
						"NAME:Color","R:=", 80,"G:=", 150,	"B:=", 250
					],
					[
						"NAME:Material",
						"Value:="		, "\"Air\""
					]
				]
			]
		])

##clen warning due to solve inside
oDesktop.ClearMessages("", "",2)

############################# Boundaries Geometry ###############################
#Circular Region Conteining Motor
oEditor.CreateCircle(
	[
		"NAME:CircleParameters","IsCovered:=", True,"XCenter:=", "0mm","YCenter:=", "0mm","ZCenter:=", "0mm","Radius:=", "%smm"%(geo["R"]),"WhichAxis:=", "Z","NumSegments:=", "0"
	], 
	[
		"NAME:Attributes","Name:=", "Region","Flags:=", "",	"Color:=", "(70 180 250)","Transparency:="	, 0.85,"PartCoordinateSystem:=","Global","UDMId:=", "","MaterialValue:=", "\"Air\"",
		"SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,	"IsMaterialEditable:=", True,"UseMaterialAppearance:=", False,"IsLightweight:="	, False
	])
#Master Boundary (lower segment)
oEditor.CreatePolyline(
	["NAME:PolylineParameters","IsPolylineCovered:=", True,"IsPolylineClosed:="	, False,
		["NAME:PolylinePoints",["NAME:PLPoint","X:=", "0mm","Y:=", "0mm","Z:=", "0mm"],
			["NAME:PLPoint","X:=", "%smm"%(geo["R"]),"Y:=", "0mm","Z:=", "0mm"]
		],
		["NAME:PolylineSegments",["NAME:PLSegment","SegmentType:=", "Line","StartIndex:=", 0,"NoOfPoints:=",2]],
		["NAME:PolylineXSection","XSectionType:=", "None","XSectionOrient:=", "Auto","XSectionWidth:=", "0mm","XSectionTopWidth:=", "0mm","XSectionHeight:=", "0mm","XSectionNumSegments:="	, "0","XSectionBendType:="	, "Corner"]
	], 
	[
		"NAME:Attributes","Name:=", "Boundary_Master","Flags:=", "","Color:=","(143 175 143)","Transparency:=", 0,"PartCoordinateSystem:=", "Global",	"UDMId:=", "","MaterialValue:=", "\"vacuum\"",
		"SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,	"IsMaterialEditable:=", True,"UseMaterialAppearance:=", False,"IsLightweight:="	, False
	])

#Slave Boundary (vertical segment)

motorangle=math.pi/geo["p"]
oEditor.CreatePolyline(
	["NAME:PolylineParameters","IsPolylineCovered:=", True,"IsPolylineClosed:="	, False,
		["NAME:PolylinePoints",["NAME:PLPoint","X:=", "0mm","Y:=", "0mm","Z:=", "0mm"],
			["NAME:PLPoint","X:=", "%smm"%(geo["R"]*math.cos(motorangle)),"Y:=", "%smm"%(geo["R"]*math.sin(motorangle)),"Z:=", "0mm"]
		],
		["NAME:PolylineSegments",["NAME:PLSegment","SegmentType:=", "Line","StartIndex:=", 0,"NoOfPoints:=",2]],
		["NAME:PolylineXSection","XSectionType:=", "None","XSectionOrient:=", "Auto","XSectionWidth:=", "0mm","XSectionTopWidth:=", "0mm","XSectionHeight:=", "0mm","XSectionNumSegments:="	, "0","XSectionBendType:="	, "Corner"]
	], 
	[
		"NAME:Attributes","Name:=", "Boundary_Slave","Flags:=", "","Color:=","(143 175 143)","Transparency:=", 0,"PartCoordinateSystem:=", "Global",	"UDMId:=", "","MaterialValue:=", "\"vacuum\"",
		"SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,	"IsMaterialEditable:=", True,"UseMaterialAppearance:=", False,"IsLightweight:="	, False
	])
#Outer Stator Arc, 0 vector potenzial
oEditor.CreatePolyline(
	["NAME:PolylineParameters","IsPolylineCovered:=",True,"IsPolylineClosed:=",False,["NAME:PolylinePoints",
			["NAME:PLPoint","X:=","%smm"%(geo["R"]),"Y:=","0mm","Z:=", "0mm"], #primo pt arco
			["NAME:PLPoint","X:=", "%smm"%(geo["R"]*math.cos(motorangle/2)),"Y:=", "%smm"%(geo["R"]*math.sin(motorangle/2)),"Z:=", "0mm"], #pt centale arco
			["NAME:PLPoint","X:=","%smm"%(geo["R"]*math.cos(motorangle)),"Y:=", "%smm"%(geo["R"]*math.sin(motorangle)),"Z:=","0mm"]#ultimo pt arco
		],
		["NAME:PolylineSegments",["NAME:PLSegment","SegmentType:=", "AngularArc","StartIndex:=", 0,	"NoOfPoints:=", 3,"NoOfSegments:=", "0",
				"ArcAngle:=", "%srad"%(motorangle),"ArcCenterX:=", "0mm","ArcCenterY:=","0mm","ArcCenterZ:=", "0mm","ArcPlane:=","XY"]],
		["NAME:PolylineXSection","XSectionType:="	, "None","XSectionOrient:="	, "Auto","XSectionWidth:=","0mm","XSectionTopWidth:=","0mm",
			"XSectionHeight:=","0mm","XSectionNumSegments:=","0","XSectionBendType:=","Corner"]
	], 
	["NAME:Attributes","Name:=","VectorPotential1","Flags:=","","Color:=","(143 175 143)","Transparency:=", 0,"PartCoordinateSystem:=", "Global",
		"UDMId:=","","MaterialValue:=", "\"vacuum\"","SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,"IsMaterialEditable:=", True,
		"UseMaterialAppearance:=", False,"IsLightweight:=", False
	])

##Region(sector) conteining moving parts (rotor)
r_mov_mid=geo["r"]+geo["g"]/2 #raggio regione mid
r_mov_out=geo["r"]+geo["g"]*3/4 #raggio regione out

#circuilar region outer rotor
oEditor.CreateCircle(
	[
		"NAME:CircleParameters","IsCovered:=", True,"XCenter:=", "0mm","YCenter:=", "0mm","ZCenter:=", "0mm","Radius:=", "%smm"%(r_mov_out),"WhichAxis:=", "Z","NumSegments:=", "0"
	], 
	[
		"NAME:Attributes","Name:=", "Rotating_band_out","Flags:=", "",	"Color:=", "(70 180 250)","Transparency:="	, 0.85,"PartCoordinateSystem:=","Global","UDMId:=", "","MaterialValue:=", "\"Air\"",
		"SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,	"IsMaterialEditable:=", True,"UseMaterialAppearance:=", False,"IsLightweight:="	, False
	])

#cut circular region
oEditor.Imprint(["NAME:Selections","Blank Parts:=","Region,Rotating_band_out","Tool Parts:=", "Boundary_Master,Boundary_Slave"	],["NAME:ImprintParameters","KeepOriginals:=",True])

facesIDs=oEditor.GetFaceIDs("Region")
faceID=int(facesIDs[1])
oEditor.DetachFaces(["NAME:Selections",	"Selections:=","Region","NewPartsModelFlag:=","Model"],["NAME:Parameters",["NAME:DetachFacesToParameters","FacesToDetach:=", [faceID]]])
oEditor.Delete(["NAME:Selections","Selections:=","Region_Detach1"])

facesIDs=oEditor.GetFaceIDs("Rotating_band_out")
faceID=int(facesIDs[1])
oEditor.DetachFaces(["NAME:Selections",	"Selections:=","Rotating_band_out","NewPartsModelFlag:=","Model"],["NAME:Parameters",["NAME:DetachFacesToParameters","FacesToDetach:=", [faceID]]])
oEditor.Delete(["NAME:Selections","Selections:=","Rotating_band_out_Detach1"])

#circuilar region mid airgap 
oEditor.CreateCircle(
	[
		"NAME:CircleParameters","IsCovered:=", True,"XCenter:=", "0mm","YCenter:=", "0mm","ZCenter:=", "0mm","Radius:=", "%smm"%(r_mov_mid),"WhichAxis:=", "Z","NumSegments:=", "0"
	], 
	[
		"NAME:Attributes","Name:=", "Rotating_band_mid","Flags:=", "",	"Color:=", "(70 180 250)","Transparency:="	, 0.85,"PartCoordinateSystem:=","Global","UDMId:=", "","MaterialValue:=", "\"Air\"",
		"SurfaceMaterialValue:=", "\"\"","SolveInside:=", True,	"IsMaterialEditable:=", True,"UseMaterialAppearance:=", False,"IsLightweight:="	, False
	])

#cut circular region
oEditor.Imprint(["NAME:Selections","Blank Parts:=","Rotating_band_mid","Tool Parts:=", "Boundary_Master,Boundary_Slave"	],["NAME:ImprintParameters","KeepOriginals:=",True])

facesIDs=oEditor.GetFaceIDs("Rotating_band_mid")
faceID=int(facesIDs[1])
oEditor.DetachFaces(["NAME:Selections",	"Selections:=","Rotating_band_mid","NewPartsModelFlag:=","Model"],["NAME:Parameters",["NAME:DetachFacesToParameters","FacesToDetach:=", [faceID]]])
oEditor.Delete(["NAME:Selections","Selections:=","Rotating_band_mid_Detach1"])


#boundaries assignment
oModule = oDesign.GetModule("BoundarySetup")
oModule.AssignMaster(["NAME:Master1","Objects:=", ["Boundary_Master"],"ReverseV:=",False])
oModule.AssignSlave(["NAME:Slave1","Objects:=",["Boundary_Slave"],"ReverseU:=", False,"Master:=","Master1","SameAsMaster:=",False])
oModule.AssignVectorPotential(["NAME:VectorPotential1","Objects:=", ["VectorPotential1"],"Value:=","0","CoordinateSystem:=", ""])

oDesktop.AddMessage ("%s"%(geo["filename"][:-4]), "Maxwell2DDesign1", 0, "end", "")

#Design Settings

oDesign.SetDesignSettings(
	[
		"NAME:Design Settings Data",
		"Perform Minimal validation:=", False,
		"EnabledObjects:="	, [],
		"PreserveTranSolnAfterDatasetEdit:=", False,
		"ComputeTransientInductance:=", False,
		"ComputeIncrementalMatrix:=", False,
		"PerfectConductorThreshold:=", 1E+30,
		"InsulatorThreshold:="	, 1,
		"ModelDepth:="		, "%smm"%(geo["l"]),
		"UseSkewModel:="	, False,
		"EnableTranTranLinkWithSimplorer:=", False,
		"BackgroundMaterialName:=", "vacuum",
		"SolveFraction:="	, False,
		"Multiplier:="		, "%d"%(geo["p"]*2)
	], 
	[
		"NAME:Model Validation Settings",
		"EntityCheckLevel:="	, "Strict",
		"IgnoreUnclassifiedObjects:=", True,
		"SkipIntersectionChecks:=", False
	])




sys.exit()
#oProject.Save()