from mbdynAdapter import MBDynHelper, MBDynAdapter
mbd = MBDynHelper()
mbd.readMsh('membrane.msh')
mbd.controlDict = {'initialTime':0,'finalTime':100,'output frequency':10}
mbd.materialDict = {'E':250, 'nu':0, 't':0.002, 'rho':500, 'C':0.000001}
adapter = MBDynAdapter(mbd)
adapter.runPreCICE()
