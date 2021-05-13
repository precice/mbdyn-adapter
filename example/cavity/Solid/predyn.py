from adapter import MBDynAdapter
from prep import MBDynPrep

mbd = MBDynPrep('membrane.msh', in_mm=False)
adapter = MBDynAdapter(mbd, debugging=False)
adapter.run_simulation()
