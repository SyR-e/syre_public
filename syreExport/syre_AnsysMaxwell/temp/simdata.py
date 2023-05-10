import pickle
import sys
simdata={
"filename":"TeslaModel3.mat", #motfilename
"filepath":"C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/", #motfilepath
"eval_type":"singt", #sim type
"iAmp":1414.000000, #Currents Amplitude
"win_avv":[[1,1,1,-3,-3,-3,2,2,2,],[1,1,1,-3,-3,-3,2,2,2,],], #matrix with windings index
"N_parallel" :1,   #
"N_cond" :6.666667e-01,   #N cond per slot
"EvalSpeed" :2000.000000,   #layer of rotor barrier
"phi_init":-1.134464, #initial current phase
"p":3,#pole pair
"q":3,#slot/(phase*p)
"radial_ribs_split":0, #area added by split barrier
"n_PM":0, #number of permanent magnet (number of magnetic segments)
"nlay" :1,   #layer of rotor barrier
"tempPP" :80.000000, #PMs temp
"gamma":0.959931,#current angle dq [rad]
"nsim":60,#simulation points
"xdeg":360.000000,#angular excursion [deg]
"corelossflag":0,#flag that activate core loss
"save_fields":0,#save fields every # steps
"dataset_counter":1,#count for dataset current
}
with open("C:/Users/paolo.ragazzo/Desktop/syre/syreExport/syre_AnsysMaxwell/temp/temp.pkl", "wb") as export:
    pickle.dump(simdata,export,protocol=2)
###
sys.exit()
