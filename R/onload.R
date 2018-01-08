'.onAttach' <- function(lib, pkg="ph2hetero")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("Loading '", desc$Package, "' version ",desc$Version);
  }
