from mbdynprecice.adapter import MBDynAdapter
from mbdynprecice.prep import MBDynPrep

mbd = MBDynPrep('membrane.msh', in_mm=False)
adapter = MBDynAdapter(mbd, debugging=False)
adapter.run_simulation()
