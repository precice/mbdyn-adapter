from mbdynadapter.adapter import MBDynAdapter
from mbdynadapter.prep import MBDynPrep

mbd = MBDynPrep('membrane.msh', in_mm=False)
adapter = MBDynAdapter(mbd, debugging=False)
adapter.run_simulation()
