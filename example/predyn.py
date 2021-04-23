from adapter import MBDynAdapter
from prep import MBDynPrep

mbd = MBDynPrep('membrane.msh')
adapter = MBDynAdapter(mbd, debugging=True, inter_mesh=True)
adapter.run_simulation()
