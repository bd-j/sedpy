import glob, os
import numpy as np
from sedpy.observate import Filter, load_filters

__all__ = ["make_md_entry"]


def make_md_entry(filternames,
                  comment="03/24: Average of normalized SCA curves",
                  directory="./prep_filters/Roman",
                  pattern="roman*.par"):


    if filternames is None:
        files = glob.glob(os.path.join(f"{directory}", f"{pattern}"))
        bases = [os.path.basename(f) for f in files]
        filternames = [b.replace(".par", "") for b in bases]

    filters = load_filters(filternames, directory=directory)
    order = np.argsort([f.wave_effective for f in filters])
    filters = np.array(filters)[order]

    for f in filters:
        name = f.name
        try:
            tel, inst, band, dummy = name.split("_")
        except(ValueError):
            tel, band = name.split("_")
            inst = ""
        weff, width = f.wave_effective, f.effective_width
        weff, width, unit = weff/10, width/10, "nm"
        row = f"|$\\textrm{{{band.upper()}}}$|{inst.upper()}|[`{name}`](filters/{name}.par)|{weff:.1f} {unit}|{width:.1f} {unit}|{comment}|"
        print(row)
