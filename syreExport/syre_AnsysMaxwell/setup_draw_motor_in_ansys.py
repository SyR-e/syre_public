import sys,os
import pickle
geodata={
"R":91.000000, #stator radius
"r":65.000000, #rotor radius
"p" :4,   #pole pair
"g" :0.800000,   #air gap thickness
"nlay" :2,   #layer of rotor barrier
"PM_CS":[[55.2387,46.1351,58.5242,54.92,],[27.5271,31.5335,19.595,10.3249,],[0.99191,0.99465,0.79113,0.63024,],[0.1269,-0.10335,0.61165,0.7764,],], #Coordinate system magnet: 0-x 1-y 2-xvect 3-yvect
"q":2,  #slot per cave per phase
"radial_ribs_split":0, #area added by split barrier
"n_PM":4, #number of permanent magnet (number of magnetic segments)
"l":156.000000, #stack length
"filepath":"C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/Eldor/", #motor"s folder
"filename":"VGB_HS.mat", #motor"s file mat name
}
#struct material names
material={
"rotor":"FeSi_3015_IronLoss",
"stator":"FeSi_3015_IronLoss",
"shaft":"ShaftAir",
"slotcond":"Copper",
"magnet":"G40EH-B",
}
with open('C:/Users/paolo.ragazzo/Desktop/syre/syreExport/syre_AnsysMaxwell/temp/temp.pkl', 'wb') as export:
    pickle.dump([geodata,material],export,protocol=2)
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64")
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64/PythonFiles/DesktopPlugin")
import ScriptEnv

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("Maxwell 2D", "Maxwell2DDesign1", "Transient", "")
oDesign = oProject.SetActiveDesign("Maxwell2DDesign1")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oProject.SaveAs("C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/Eldor/VGB_HS.aedt", True)
oDesktop.ClearMessages("", "",3)
#########Add materials to library
#########Air##########
oDefinitionManager = oProject.GetDefinitionManager()
oDefinitionManager.AddMaterial(
	[
		"NAME:Air ",
		"CoordinateSystemType:=", "Cartesian",
		"BulkOrSurfaceType:="	, 1,
		[
			"NAME:PhysicsTypes",
			"set:="			, ["Electromagnetic"]
		]
	])
############Magnet############
oDefinitionManager.AddMaterial(
	[
	"NAME:G40EH-B",
	"CoordinateSystemType:=", "Cartesian",
	"BulkOrSurfaceType:="	, 1,
	[
		"NAME:PhysicsTypes",
		"set:="			, ["Electromagnetic"]
	],
	[
		"NAME:ModifierData",
		[
			"NAME:ThermalModifierData",
			"modifier_data:="	, "thermal_modifier_data",
			[
				"NAME:all_thermal_modifiers",
				[
					"NAME:one_thermal_modifier",
					"Property::="		, "magnetic_coercivity",
					"Index::="		, 0,
					"prop_modifier:="	, "thermal_modifier",
					"use_free_form:="	, False,
					"Tref:="		, "20cel",
					"C1:="			, "-0.001200",
					"C2:="			, "0",
					"TL:="			, "-273.15cel",
					"TU:="			, "1000cel",
					"auto_calculation:="	, True
				]
			]
		]
	],
	"permeability:="	, "1.050000",
	"conductivity:="	, "714290.000000",
	[
		"NAME:magnetic_coercivity",
		"property_type:="	, "VectorProperty",
		"Magnitude:="		, "-907183.175624A_per_meter",
		"DirComp1:="		, "1",
		"DirComp2:="		, "0",
		"DirComp3:="		, "0"
	],
	"mass_density:="	, "7600.000000"
])
########## Slot conductor ##########
oDefinitionManager.AddMaterial(
[
	"NAME:Copper",
	"CoordinateSystemType:=", "Cartesian",
	"BulkOrSurfaceType:="	, 1,
	[
		"NAME:PhysicsTypes",
		"set:="			, ["Electromagnetic","Thermal"]
	],
	"conductivity:="	, "58000000.000000",
	"thermal_conductivity:=", "400.000000",
	"mass_density:="	, "8940.000000"
])
#########Add materials to library
#########ShaftAir##########
oDefinitionManager = oProject.GetDefinitionManager()
oDefinitionManager.AddMaterial(
	[
		"NAME:ShaftAir ",
		"CoordinateSystemType:=", "Cartesian",
		"BulkOrSurfaceType:="	, 1,
		[
			"NAME:PhysicsTypes",
			"set:="			, ["Electromagnetic"]
		]
	])
############ Rotor"s plates ###########
oDefinitionManager.AddMaterial(
[
	"NAME:FeSi_3015_IronLoss",
	"CoordinateSystemType:=", "Cartesian",
	"BulkOrSurfaceType:="	, 1,
	[
		"NAME:PhysicsTypes",
		"set:="			, ["Electromagnetic"]
	],
	[
		"NAME:permeability",
		"property_type:="	, "nonlinear",
		"BTypeForSingleCurve:="	, "normal",
		"HUnit:="		, "A_per_meter",
		"BUnit:="		, "tesla",
		"IsTemperatureDependent:=", False,
		[
			"NAME:BHCoordinates",
			[
				"NAME:DimUnits", 
				"", 
				""
			],
			["NAME:Coordinate",["NAME:CoordPoint",0.000000,0.000000]],
			["NAME:Coordinate",["NAME:CoordPoint",10.000000,0.008000]],
			["NAME:Coordinate",["NAME:CoordPoint",20.000000,0.045000]],
			["NAME:Coordinate",["NAME:CoordPoint",30.000000,0.101000]],
			["NAME:Coordinate",["NAME:CoordPoint",40.000000,0.231100]],
			["NAME:Coordinate",["NAME:CoordPoint",50.000000,0.439100]],
			["NAME:Coordinate",["NAME:CoordPoint",60.000000,0.637100]],
			["NAME:Coordinate",["NAME:CoordPoint",70.000000,0.780100]],
			["NAME:Coordinate",["NAME:CoordPoint",80.000000,0.882100]],
			["NAME:Coordinate",["NAME:CoordPoint",90.000000,0.961100]],
			["NAME:Coordinate",["NAME:CoordPoint",100.000000,1.020100]],
			["NAME:Coordinate",["NAME:CoordPoint",125.000000,1.119200]],
			["NAME:Coordinate",["NAME:CoordPoint",150.000000,1.180200]],
			["NAME:Coordinate",["NAME:CoordPoint",200.000000,1.252300]],
			["NAME:Coordinate",["NAME:CoordPoint",250.000000,1.293300]],
			["NAME:Coordinate",["NAME:CoordPoint",350.000000,1.342400]],
			["NAME:Coordinate",["NAME:CoordPoint",500.000000,1.381600]],
			["NAME:Coordinate",["NAME:CoordPoint",750.000000,1.418900]],
			["NAME:Coordinate",["NAME:CoordPoint",1000.000000,1.445300]],
			["NAME:Coordinate",["NAME:CoordPoint",1250.000000,1.464600]],
			["NAME:Coordinate",["NAME:CoordPoint",1500.000000,1.481900]],
			["NAME:Coordinate",["NAME:CoordPoint",2000.000000,1.511500]],
			["NAME:Coordinate",["NAME:CoordPoint",2500.000000,1.536100]],
			["NAME:Coordinate",["NAME:CoordPoint",5000.000000,1.631300]],
			["NAME:Coordinate",["NAME:CoordPoint",7500.000000,1.706400]],
			["NAME:Coordinate",["NAME:CoordPoint",10000.000000,1.763600]],
			["NAME:Coordinate",["NAME:CoordPoint",12000.000000,1.800000]],
			["NAME:Coordinate",["NAME:CoordPoint",15000.000000,1.838900]],
			["NAME:Coordinate",["NAME:CoordPoint",20000.000000,1.881800]],
			["NAME:Coordinate",["NAME:CoordPoint",30000.000000,1.932500]],
			["NAME:Coordinate",["NAME:CoordPoint",40000.000000,1.964800]],
			["NAME:Coordinate",["NAME:CoordPoint",50000.000000,1.989300]],
			["NAME:Coordinate",["NAME:CoordPoint",60000.000000,2.010000]],
			["NAME:Coordinate",["NAME:CoordPoint",70000.000000,2.028300]],
			["NAME:Coordinate",["NAME:CoordPoint",80000.000000,2.045300]],
			["NAME:Coordinate",["NAME:CoordPoint",90000.000000,2.061300]],
			["NAME:Coordinate",["NAME:CoordPoint",100000.000000,2.076600]],
			["NAME:Coordinate",["NAME:CoordPoint",200000.000000,2.214700]],
			["NAME:Coordinate",["NAME:CoordPoint",300000.000000,2.344600]],
			["NAME:Coordinate",["NAME:CoordPoint",400000.000000,2.472300]],
			["NAME:Coordinate",["NAME:CoordPoint",500000.000000,2.599300]],
			["NAME:Coordinate",["NAME:CoordPoint",600000.000000,2.725800]],
			["NAME:Coordinate",["NAME:CoordPoint",700000.000000,2.852000]],
			["NAME:Coordinate",["NAME:CoordPoint",800000.000000,2.978100]],
			["NAME:Coordinate",["NAME:CoordPoint",900000.000000,3.104200]],
			["NAME:Coordinate",["NAME:CoordPoint",1000000.000000,3.230100]],
		],
		[
			"NAME:Temperatures"
		]
	],
	[
		"NAME:magnetic_coercivity",
		"property_type:="	, "VectorProperty",
		"Magnitude:="		, "0A_per_meter",
		"DirComp1:="		, "1",
		"DirComp2:="		, "0",
		"DirComp3:="		, "0"
	],
	[
		"NAME:core_loss_type",
		"property_type:="	, "ChoiceProperty",
		"Choice:="		, "Electrical Steel"
	],
	"core_loss_kh:="	, "135.027726",
	"core_loss_kc:="	, "0.044536",
	"core_loss_ke:="	, "7.750784",
	"core_loss_kdc:="	, "0",
	"mass_density:="	, "7600.000000",
	"core_loss_equiv_cut_depth:=", "0.001meter"
])

sys.exit()
