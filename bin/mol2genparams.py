#!/usr/bin/env python

import sys
import vscreenml_v2

if __name__ == "__main__":

    genpot_path = vscreenml_v2.__file__
    genpot_path = genpot_path[0:genpot_path.rfind("/")]
    genpot_path += "/generic_potential/"

    sys.path.append(genpot_path)

    from vscreenml_v2.generic_potential.mol2genparams import main
    from vscreenml_v2.generic_potential.BasicClasses import OptionClass

    option = OptionClass(sys.argv)
    main(option)