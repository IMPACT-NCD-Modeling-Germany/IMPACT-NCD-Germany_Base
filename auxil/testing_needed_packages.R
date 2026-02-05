needed_packages <- yaml::read_yaml("./dependencies.yaml")
installed_packages <- installed.packages()[, "Package"]
already_installed_packages <- intersect(needed_packages, installed_packages)
missing_packages <- setdiff(needed_packages, already_installed_packages)
