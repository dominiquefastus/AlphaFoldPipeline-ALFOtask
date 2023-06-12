#!/usr/bin/env coot
# python script for coot - generated by dimple
set_nomenclature_errors_on_read("ignore")
import inspect, os
this_file = inspect.getfile(inspect.currentframe())
this_dir = os.path.dirname(os.path.abspath(this_file))
molecule = read_pdb(os.path.join(this_dir, "final.pdb"))
set_rotation_centre(-34.21, 16.3, -11.19)
set_zoom(30.)
set_view_quaternion(-0.260911, -0.600234, 0, 0.756071)
refl = os.path.join(this_dir, "final.mtz")
map21 = make_and_draw_map(refl, "FWT", "PHWT", "", 0, 0)
map11 = make_and_draw_map(refl, "DELFWT", "PHDELWT", "", 0, 1)