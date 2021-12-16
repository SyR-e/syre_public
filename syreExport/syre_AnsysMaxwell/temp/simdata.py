import pickle
import sys
simdata={
"filename":"VGB_LS.mat", #motfilename
"filepath":"C:/Users/paolo.ragazzo/Desktop/syre/motorExamples/Eldor/", #motfilepath
"eval_type":"singtIron", #sim type
"iAmp":500.020000, #Currents Amplitude
"win_avv":[[1,1,-3,-3,2,2,],[1,1,-3,-3,2,2,],], #matrix with windings index
"N_parallel" :1,   #
"N_cond" :1,   #N cond per slot
"EvalSpeed" :13000.000000,   #layer of rotor barrier
"phi_init":-1.315629, #initial current phase
"p":4,#pole pair
"q":2,#slot/(phase*p)
"radial_ribs_split":0, #area added by split barrier
"n_PM":0, #number of permanent magnet (number of magnetic segments)
"nlay" :2,   #layer of rotor barrier
"tempPP" :80.000000, #PMs temp
"gamma":0.778766,#current angle dq [rad]
"nsim":50,#simulation points
"xdeg":420.000000,#angular excursion [deg]
"corelossflag":1,#flag that activate core loss
"save_fields":0,#save fields every # steps
"dataset_counter":1,#count for dataset current
}
with open("C:/Users/paolo.ragazzo/Desktop/syre/syreExport/syre_AnsysMaxwell/temp/temp.pkl", "wb") as export:
    pickle.dump(simdata,export,protocol=2)
###
sys.exit()
