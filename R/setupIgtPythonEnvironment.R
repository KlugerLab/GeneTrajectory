#' Set up Python environment for running IGT (Inference of Gene Trajectories) using reticulate package in R
#'
#' @export
#'

setupIgtPythonEnvironment <- function(condaenv='igt',
                                      virtualenv='igt',
                                      packages=c('numpy', 'scipy', 'pot', 'tqdm'),
                                      create=FALSE){
  if(reticulate::py_available()){
    message(paste('Already set up to ', reticulate::py_exe()))
    return(TRUE)
  }
  
  if(reticulate::condaenv_exists(condaenv)){
    message(paste('Using conda environment', condaenv))
    reticulate::use_condaenv(condaenv = condaenv)
    return(TRUE)
  }
  
  if(reticulate::virtualenv_exists(virtualenv)){
    message(paste('Using virtualenv environment', condaenv))
    reticulate::use_virtualenv(virtualenv=virtualenv)
    return(TRUE)
  }
  
  
  if(create & ! reticulate::condaenv_exists(condaenv)){
    reticulate::conda_create(condaenv, packages = packages)
    reticulate::use_condaenv(condaenv = condaenv)
    return(TRUE)
  }
  
  if(create & ! reticulate::virtualenv_exists(virtualenv)){
    virtualenv_create(virtualenv, packages=packages)
    reticulate::use_virtualenv(virtualenv=virtualenv)
    return(TRUE)
  }
  
  return(FALSE)
}


