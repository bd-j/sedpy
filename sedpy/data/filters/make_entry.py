import glob, os
import numpy as np
from sedpy.observate import Filter, load_filters

dirn = "./prep_filters/Roman"
files = glob.glob(f"{dirn}/*roman*par")
bases = [os.path.basename(f) for f in files]
names = [b.replace(".par", "") for b in bases]

filters = load_filters(names, directory=dirn)
order = np.argsort([f.wave_effective for f in filters])
filters = filters[order]

comment = "03/24: Average of normalized SCA curves"

for f in filters:
    name = f.name
    tel, inst, band = name.split("_")
    weff, width = f.wave_effective, f.effective_width
    weff, width, unit = weff/10, width/10, "nm"
    row = f"|$\\textrm{{{band.upper()}}}$|{inst.upper()}|[`{name}`](filters/{name}.par)|{weff:.1f} {unit}|{width:.1f} {unit}|{comment}|"
    print(row)