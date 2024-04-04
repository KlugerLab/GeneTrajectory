#' Check the Python environment for running IGT (Inference of Gene Trajectories) using reticulate package in R
#'
#' @export
#'

checkIgtPythonEnvironment <- function(python.packages = c('ot', 'tqdm', 'scipy.io'), check.numpy=TRUE){
  if(!reticulate::py_available()){
    setupIgtPythonEnvironment()
    reticulate::py_available(initialize=TRUE)
  }
  if(!reticulate::py_available()){
    stop('Python not available. Please install Python, e.g. by using using reticulate::install_python()')
  }
  if(check.numpy & ! reticulate::py_numpy_available()){
    stop(sprintf('Package "numpy" not available.'))
  }
  for (p in python.packages) {
    if(! reticulate::py_module_available(p)){
      stop(sprintf('Package "%s" not available.', p))
    }
  }
}