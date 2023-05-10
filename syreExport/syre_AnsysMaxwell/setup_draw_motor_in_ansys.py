import sys,os
import pickle
geodata={
"R":112.500000, #stator radius
"r":74.950000, #rotor radius
"p" :3,   #pole pair
"g" :0.700000,   #air gap thickness
"nlay" :1,   #layer of rotor barrier
"PM_CS":[[46.0023,59.9243,],[42.6352,18.5216,],[0.98942,0.62033,],[0.14506,0.78434,],], #Coordinate system magnet: 0-x 1-y 2-xvect 3-yvect
"q":3,  #slot per cave per phase
"radial_ribs_split":0, #area added by split barrier
"n_PM":2, #number of permanent magnet (number of magnetic segments)
"l":134.000000, #stack length
"filepath":"C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/", #motor"s folder
"filename":"TeslaModel3.mat", #motor"s file mat name
}
#struct material names
material={
"rotor":"M270-35A",
"stator":"M270-35A",
"shaft":"ShaftAir",
"slotcond":"Copper",
"magnet":"BMN-52UH",
}
with open('C:/Users/paolo.ragazzo/Desktop/syre/syreExport/syre_AnsysMaxwell/temp/temp.pkl', 'wb') as export:
    pickle.dump([geodata,material],export,protocol=2)
sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64")
sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64/PythonFiles/DesktopPlugin")
import ScriptEnv

ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oProject = oDesktop.NewProject()
oProject.InsertDesign("Maxwell 2D", "Maxwell2DDesign1", "Transient", "")
oDesign = oProject.SetActiveDesign("Maxwell2DDesign1")
oEditor = oDesign.SetActiveEditor("3D Modeler")
oProject.SaveAs("C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/TeslaModel3.aedt", True)
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
	"NAME:BMN-52UH",
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
					"C1:="			, "-0.001150",
					"C2:="			, "0",
					"TL:="			, "-273.15cel",
					"TU:="			, "1000cel",
					"auto_calculation:="	, True
				]
			]
		]
	],
	"permeability:="	, "1.050000",
	"conductivity:="	, "666000.000000",
	[
		"NAME:magnetic_coercivity",
		"property_type:="	, "VectorProperty",
		"Magnitude:="		, "-1098926.988015A_per_meter",
		"DirComp1:="		, "1",
		"DirComp2:="		, "0",
		"DirComp3:="		, "0"
	],
	"mass_density:="	, "7550.000000"
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
	"NAME:M270-35A",
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
			["NAME:Coordinate",["NAME:CoordPoint",3.276787,0.004542]],
			["NAME:Coordinate",["NAME:CoordPoint",6.114593,0.009035]],
			["NAME:Coordinate",["NAME:CoordPoint",8.460242,0.013222]],
			["NAME:Coordinate",["NAME:CoordPoint",10.462328,0.017196]],
			["NAME:Coordinate",["NAME:CoordPoint",12.199267,0.020991]],
			["NAME:Coordinate",["NAME:CoordPoint",13.729824,0.024640]],
			["NAME:Coordinate",["NAME:CoordPoint",15.095619,0.028169]],
			["NAME:Coordinate",["NAME:CoordPoint",16.327217,0.031597]],
			["NAME:Coordinate",["NAME:CoordPoint",17.447600,0.034941]],
			["NAME:Coordinate",["NAME:CoordPoint",18.474440,0.038213]],
			["NAME:Coordinate",["NAME:CoordPoint",19.421602,0.041423]],
			["NAME:Coordinate",["NAME:CoordPoint",20.300151,0.044579]],
			["NAME:Coordinate",["NAME:CoordPoint",21.119056,0.047689]],
			["NAME:Coordinate",["NAME:CoordPoint",21.885681,0.050758]],
			["NAME:Coordinate",["NAME:CoordPoint",22.606148,0.053792]],
			["NAME:Coordinate",["NAME:CoordPoint",23.285597,0.056795]],
			["NAME:Coordinate",["NAME:CoordPoint",23.928387,0.059771]],
			["NAME:Coordinate",["NAME:CoordPoint",24.538245,0.062724]],
			["NAME:Coordinate",["NAME:CoordPoint",25.118382,0.065655]],
			["NAME:Coordinate",["NAME:CoordPoint",25.671584,0.068569]],
			["NAME:Coordinate",["NAME:CoordPoint",26.200284,0.071468]],
			["NAME:Coordinate",["NAME:CoordPoint",26.706619,0.074353]],
			["NAME:Coordinate",["NAME:CoordPoint",27.192476,0.077227]],
			["NAME:Coordinate",["NAME:CoordPoint",27.659531,0.080091]],
			["NAME:Coordinate",["NAME:CoordPoint",28.109277,0.082948]],
			["NAME:Coordinate",["NAME:CoordPoint",28.543052,0.085800]],
			["NAME:Coordinate",["NAME:CoordPoint",28.962059,0.088646]],
			["NAME:Coordinate",["NAME:CoordPoint",29.367384,0.091490]],
			["NAME:Coordinate",["NAME:CoordPoint",29.760012,0.094332]],
			["NAME:Coordinate",["NAME:CoordPoint",30.140836,0.097174]],
			["NAME:Coordinate",["NAME:CoordPoint",30.510673,0.100016]],
			["NAME:Coordinate",["NAME:CoordPoint",30.870270,0.102861]],
			["NAME:Coordinate",["NAME:CoordPoint",31.220314,0.105709]],
			["NAME:Coordinate",["NAME:CoordPoint",31.561434,0.108561]],
			["NAME:Coordinate",["NAME:CoordPoint",31.894215,0.111419]],
			["NAME:Coordinate",["NAME:CoordPoint",32.219196,0.114283]],
			["NAME:Coordinate",["NAME:CoordPoint",32.536879,0.117154]],
			["NAME:Coordinate",["NAME:CoordPoint",32.847730,0.120035]],
			["NAME:Coordinate",["NAME:CoordPoint",33.152185,0.122924]],
			["NAME:Coordinate",["NAME:CoordPoint",33.450652,0.125825]],
			["NAME:Coordinate",["NAME:CoordPoint",33.743514,0.128737]],
			["NAME:Coordinate",["NAME:CoordPoint",34.031132,0.131661]],
			["NAME:Coordinate",["NAME:CoordPoint",34.313846,0.134600]],
			["NAME:Coordinate",["NAME:CoordPoint",34.591978,0.137553]],
			["NAME:Coordinate",["NAME:CoordPoint",34.865834,0.140522]],
			["NAME:Coordinate",["NAME:CoordPoint",35.135705,0.143507]],
			["NAME:Coordinate",["NAME:CoordPoint",35.401869,0.146511]],
			["NAME:Coordinate",["NAME:CoordPoint",35.664594,0.149533]],
			["NAME:Coordinate",["NAME:CoordPoint",35.924134,0.152576]],
			["NAME:Coordinate",["NAME:CoordPoint",36.180736,0.155640]],
			["NAME:Coordinate",["NAME:CoordPoint",36.434640,0.158726]],
			["NAME:Coordinate",["NAME:CoordPoint",36.686077,0.161836]],
			["NAME:Coordinate",["NAME:CoordPoint",36.935274,0.164972]],
			["NAME:Coordinate",["NAME:CoordPoint",37.182450,0.168133]],
			["NAME:Coordinate",["NAME:CoordPoint",37.427824,0.171322]],
			["NAME:Coordinate",["NAME:CoordPoint",37.671610,0.174540]],
			["NAME:Coordinate",["NAME:CoordPoint",37.914018,0.177789]],
			["NAME:Coordinate",["NAME:CoordPoint",38.155260,0.181070]],
			["NAME:Coordinate",["NAME:CoordPoint",38.395547,0.184384]],
			["NAME:Coordinate",["NAME:CoordPoint",38.635088,0.187734]],
			["NAME:Coordinate",["NAME:CoordPoint",38.874097,0.191121]],
			["NAME:Coordinate",["NAME:CoordPoint",39.112790,0.194547]],
			["NAME:Coordinate",["NAME:CoordPoint",39.351384,0.198014]],
			["NAME:Coordinate",["NAME:CoordPoint",39.590104,0.201524]],
			["NAME:Coordinate",["NAME:CoordPoint",39.829179,0.205079]],
			["NAME:Coordinate",["NAME:CoordPoint",40.068846,0.208681]],
			["NAME:Coordinate",["NAME:CoordPoint",40.309351,0.212334]],
			["NAME:Coordinate",["NAME:CoordPoint",40.550950,0.216039]],
			["NAME:Coordinate",["NAME:CoordPoint",40.793910,0.219799]],
			["NAME:Coordinate",["NAME:CoordPoint",41.038513,0.223617]],
			["NAME:Coordinate",["NAME:CoordPoint",41.285056,0.227497]],
			["NAME:Coordinate",["NAME:CoordPoint",41.533854,0.231442]],
			["NAME:Coordinate",["NAME:CoordPoint",41.785244,0.235454]],
			["NAME:Coordinate",["NAME:CoordPoint",42.039585,0.239539]],
			["NAME:Coordinate",["NAME:CoordPoint",42.297265,0.243699]],
			["NAME:Coordinate",["NAME:CoordPoint",42.558702,0.247941]],
			["NAME:Coordinate",["NAME:CoordPoint",42.824350,0.252267]],
			["NAME:Coordinate",["NAME:CoordPoint",43.094702,0.256684]],
			["NAME:Coordinate",["NAME:CoordPoint",43.370301,0.261196]],
			["NAME:Coordinate",["NAME:CoordPoint",43.651739,0.265811]],
			["NAME:Coordinate",["NAME:CoordPoint",43.940203,0.270542]],
			["NAME:Coordinate",["NAME:CoordPoint",44.238730,0.275438]],
			["NAME:Coordinate",["NAME:CoordPoint",44.548036,0.280510]],
			["NAME:Coordinate",["NAME:CoordPoint",44.868763,0.285768]],
			["NAME:Coordinate",["NAME:CoordPoint",45.201605,0.291223]],
			["NAME:Coordinate",["NAME:CoordPoint",45.547319,0.296887]],
			["NAME:Coordinate",["NAME:CoordPoint",45.906732,0.302770]],
			["NAME:Coordinate",["NAME:CoordPoint",46.280750,0.308889]],
			["NAME:Coordinate",["NAME:CoordPoint",46.670366,0.315257]],
			["NAME:Coordinate",["NAME:CoordPoint",47.076677,0.321890]],
			["NAME:Coordinate",["NAME:CoordPoint",47.500893,0.328806]],
			["NAME:Coordinate",["NAME:CoordPoint",47.944354,0.336025]],
			["NAME:Coordinate",["NAME:CoordPoint",48.408551,0.343567]],
			["NAME:Coordinate",["NAME:CoordPoint",48.895149,0.351457]],
			["NAME:Coordinate",["NAME:CoordPoint",49.406012,0.359720]],
			["NAME:Coordinate",["NAME:CoordPoint",49.943240,0.368385]],
			["NAME:Coordinate",["NAME:CoordPoint",50.509204,0.377485]],
			["NAME:Coordinate",["NAME:CoordPoint",51.106596,0.387054]],
			["NAME:Coordinate",["NAME:CoordPoint",51.738489,0.397133]],
			["NAME:Coordinate",["NAME:CoordPoint",52.408410,0.407766]],
			["NAME:Coordinate",["NAME:CoordPoint",53.120425,0.419005]],
			["NAME:Coordinate",["NAME:CoordPoint",53.879260,0.430908]],
			["NAME:Coordinate",["NAME:CoordPoint",54.690431,0.443538]],
			["NAME:Coordinate",["NAME:CoordPoint",55.560433,0.456972]],
			["NAME:Coordinate",["NAME:CoordPoint",56.496959,0.471295]],
			["NAME:Coordinate",["NAME:CoordPoint",57.509196,0.486604]],
			["NAME:Coordinate",["NAME:CoordPoint",58.608200,0.503015]],
			["NAME:Coordinate",["NAME:CoordPoint",59.807376,0.520657]],
			["NAME:Coordinate",["NAME:CoordPoint",61.123111,0.539683]],
			["NAME:Coordinate",["NAME:CoordPoint",62.575568,0.560265]],
			["NAME:Coordinate",["NAME:CoordPoint",64.189699,0.582603]],
			["NAME:Coordinate",["NAME:CoordPoint",65.996464,0.606916]],
			["NAME:Coordinate",["NAME:CoordPoint",68.034202,0.633441]],
			["NAME:Coordinate",["NAME:CoordPoint",70.349880,0.662415]],
			["NAME:Coordinate",["NAME:CoordPoint",72.999623,0.694037]],
			["NAME:Coordinate",["NAME:CoordPoint",76.047205,0.728403]],
			["NAME:Coordinate",["NAME:CoordPoint",79.558510,0.765409]],
			["NAME:Coordinate",["NAME:CoordPoint",83.589963,0.804631]],
			["NAME:Coordinate",["NAME:CoordPoint",88.171938,0.845245]],
			["NAME:Coordinate",["NAME:CoordPoint",93.294735,0.886100]],
			["NAME:Coordinate",["NAME:CoordPoint",98.908944,0.925950]],
			["NAME:Coordinate",["NAME:CoordPoint",104.943590,0.963771]],
			["NAME:Coordinate",["NAME:CoordPoint",111.330056,0.998940]],
			["NAME:Coordinate",["NAME:CoordPoint",118.017431,1.031222]],
			["NAME:Coordinate",["NAME:CoordPoint",124.976096,1.060658]],
			["NAME:Coordinate",["NAME:CoordPoint",132.194547,1.087441]],
			["NAME:Coordinate",["NAME:CoordPoint",139.674490,1.111822]],
			["NAME:Coordinate",["NAME:CoordPoint",147.426649,1.134065]],
			["NAME:Coordinate",["NAME:CoordPoint",155.467825,1.154418]],
			["NAME:Coordinate",["NAME:CoordPoint",163.819053,1.173103]],
			["NAME:Coordinate",["NAME:CoordPoint",172.504540,1.190318]],
			["NAME:Coordinate",["NAME:CoordPoint",181.551128,1.206233]],
			["NAME:Coordinate",["NAME:CoordPoint",190.988105,1.220994]],
			["NAME:Coordinate",["NAME:CoordPoint",200.847236,1.234731]],
			["NAME:Coordinate",["NAME:CoordPoint",211.162964,1.247553]],
			["NAME:Coordinate",["NAME:CoordPoint",221.972737,1.259557]],
			["NAME:Coordinate",["NAME:CoordPoint",233.317449,1.270827]],
			["NAME:Coordinate",["NAME:CoordPoint",245.241993,1.281438]],
			["NAME:Coordinate",["NAME:CoordPoint",257.795946,1.291454]],
			["NAME:Coordinate",["NAME:CoordPoint",271.034398,1.300935]],
			["NAME:Coordinate",["NAME:CoordPoint",285.018964,1.309934]],
			["NAME:Coordinate",["NAME:CoordPoint",299.819023,1.318498]],
			["NAME:Coordinate",["NAME:CoordPoint",315.513246,1.326671]],
			["NAME:Coordinate",["NAME:CoordPoint",332.191496,1.334492]],
			["NAME:Coordinate",["NAME:CoordPoint",349.957213,1.342001]],
			["NAME:Coordinate",["NAME:CoordPoint",368.930446,1.349231]],
			["NAME:Coordinate",["NAME:CoordPoint",389.251734,1.356218]],
			["NAME:Coordinate",["NAME:CoordPoint",411.087165,1.362993]],
			["NAME:Coordinate",["NAME:CoordPoint",434.635037,1.369591]],
			["NAME:Coordinate",["NAME:CoordPoint",460.134770,1.376045]],
			["NAME:Coordinate",["NAME:CoordPoint",487.879039,1.382389]],
			["NAME:Coordinate",["NAME:CoordPoint",518.230610,1.388661]],
			["NAME:Coordinate",["NAME:CoordPoint",551.646186,1.394903]],
			["NAME:Coordinate",["NAME:CoordPoint",588.711043,1.401161]],
			["NAME:Coordinate",["NAME:CoordPoint",630.190743,1.407491]],
			["NAME:Coordinate",["NAME:CoordPoint",677.110879,1.413963]],
			["NAME:Coordinate",["NAME:CoordPoint",730.884830,1.420665]],
			["NAME:Coordinate",["NAME:CoordPoint",793.527840,1.427717]],
			["NAME:Coordinate",["NAME:CoordPoint",868.035803,1.435288]],
			["NAME:Coordinate",["NAME:CoordPoint",959.101751,1.443633]],
			["NAME:Coordinate",["NAME:CoordPoint",1074.589186,1.453164]],
			["NAME:Coordinate",["NAME:CoordPoint",1228.903408,1.464609]],
			["NAME:Coordinate",["NAME:CoordPoint",1451.854843,1.479414]],
			["NAME:Coordinate",["NAME:CoordPoint",1812.597843,1.500730]],
			["NAME:Coordinate",["NAME:CoordPoint",2396.292278,1.531321]],
			["NAME:Coordinate",["NAME:CoordPoint",3340.729712,1.574972]],
			["NAME:Coordinate",["NAME:CoordPoint",4868.861582,1.635173]],
			["NAME:Coordinate",["NAME:CoordPoint",6056.786073,1.674606]],
			["NAME:Coordinate",["NAME:CoordPoint",7126.730120,1.705354]],
			["NAME:Coordinate",["NAME:CoordPoint",8135.946630,1.730750]],
			["NAME:Coordinate",["NAME:CoordPoint",9119.387574,1.752580]],
			["NAME:Coordinate",["NAME:CoordPoint",10098.258480,1.771845]],
			["NAME:Coordinate",["NAME:CoordPoint",11087.221330,1.789166]],
			["NAME:Coordinate",["NAME:CoordPoint",12097.589930,1.804957]],
			["NAME:Coordinate",["NAME:CoordPoint",13138.977850,1.819514]],
			["NAME:Coordinate",["NAME:CoordPoint",14220.250510,1.833057]],
			["NAME:Coordinate",["NAME:CoordPoint",15350.151270,1.845755]],
			["NAME:Coordinate",["NAME:CoordPoint",16537.778850,1.857747]],
			["NAME:Coordinate",["NAME:CoordPoint",17793.012470,1.869145]],
			["NAME:Coordinate",["NAME:CoordPoint",19126.948680,1.880049]],
			["NAME:Coordinate",["NAME:CoordPoint",20552.403480,1.890545]],
			["NAME:Coordinate",["NAME:CoordPoint",22084.537740,1.900713]],
			["NAME:Coordinate",["NAME:CoordPoint",23741.681730,1.910631]],
			["NAME:Coordinate",["NAME:CoordPoint",25546.469320,1.920376]],
			["NAME:Coordinate",["NAME:CoordPoint",27527.455720,1.930028]],
			["NAME:Coordinate",["NAME:CoordPoint",29721.505960,1.939676]],
			["NAME:Coordinate",["NAME:CoordPoint",32177.452920,1.949423]],
			["NAME:Coordinate",["NAME:CoordPoint",34961.935420,1.959397]],
			["NAME:Coordinate",["NAME:CoordPoint",38169.180130,1.969762]],
			["NAME:Coordinate",["NAME:CoordPoint",41938.394230,1.980747]],
			["NAME:Coordinate",["NAME:CoordPoint",46487.096960,1.992692]],
			["NAME:Coordinate",["NAME:CoordPoint",52181.595710,2.006150]],
			["NAME:Coordinate",["NAME:CoordPoint",59707.652570,2.022124]],
			["NAME:Coordinate",["NAME:CoordPoint",70576.241300,2.042772]],
			["NAME:Coordinate",["NAME:CoordPoint",88161.987270,2.072563]],
			["NAME:Coordinate",["NAME:CoordPoint",116616.322000,2.115877]],
			["NAME:Coordinate",["NAME:CoordPoint",162656.402600,2.180377]],
			["NAME:Coordinate",["NAME:CoordPoint",237150.818000,2.279283]],
			["NAME:Coordinate",["NAME:CoordPoint",310949.486200,2.374766]],
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
	"core_loss_kh:="	, "167.239114",
	"core_loss_kc:="	, "0.387508",
	"core_loss_ke:="	, "0.034983",
	"core_loss_kdc:="	, "0",
	"mass_density:="	, "7650.000000",
	"core_loss_equiv_cut_depth:=", "0.001meter"
])

sys.exit()
