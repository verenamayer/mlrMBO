# magic mboTemplate - in this function the mbo magic for all our mbo approaches
# does happen - model fitting and point proposal in a generall way. the respective
# mbo algorithms differ in the subfunctions.
# - usually the mboTemplate is started from an OptProblem which is an environment.
# - mboTemplate can also be called from mboContinue from a saved OptState
# - The opt.state is also en environment linking to the main Objects
#    - OptProblem (constant; stores the information which define the problem)
#    - OptPath (stores all information about function evaluations)
#    - OptResult (stores information which should be part of the later constructed mboResult)
#    - (see respective source files for further information)

mboTemplate = function(obj) {
  UseMethod("mboTemplate")
}

# Creates the initial OptState and then runs the template on it
mboTemplate.OptProblem = function(obj) {
  opt.state = makeOptState(obj)
  # evaluate initial design (if y not given) and log to optpath
  evalMBODesign.OptState(opt.state)
  finalizeMboLoop(opt.state)
  
  # if you want epsilon distance for proposed points, set epsilon start value
  if (opt.state$opt.problem$control$infill.eps.proposed.points.rf) {
    if (opt.state$opt.problem$control$infill.eps.dist.proposed.points == "gower") {
      dist_designpoints = gower.dist(data.x = opt.state$opt.problem$design)
      # without diag points (distance of the point itself)
      diag(dist_designpoints) = NA
      dist_designpoints = dist_designpoints[(which(!is.na(dist_designpoints)))]
    } else if (opt.state$opt.problem$control$infill.eps.dist.proposed.points == "euclidean") {
      dist_designpoints = dist(opt.state$opt.problem$design, method = "euclidean")
    } else {
      stop("You can only set gower or euclidean for eps.dist.proposed.points!")
    }
    # For now, set default epsilon start here; for example: mean distance designpoints / 2
    if (opt.state$opt.problem$control$infill.eps.start.proposed.points.rf == "mean.dist") {
      opt.state$opt.problem$control$infill.eps.start.proposed.points.rf = mean(dist_designpoints) / 2
    }
    if (is.numeric(opt.state$opt.problem$control$infill.eps.start.proposed.points.rf) != TRUE) {
      stop("eps_start must be numeric!")
    }
  }
  
  mboTemplate(opt.state)
}

# Runs the mbo iterations on any given OptState until termination criterion is fulfilled
mboTemplate.OptState = function(obj) {
  opt.state = obj
  setOptStateLoopStarttime(opt.state)
  # check if budget is already exceeded after intitial design creation
  terminate = getOptStateTermination(opt.state)
  if (terminate$term) {
    opt.problem = getOptStateOptProblem(opt.state)
    showInfo(getOptProblemShowInfo(opt.problem), "%s. The stopping conditions
      was satisfied right after the creation of the initial design!", terminate$message)
    return(opt.state)
  }
  
  # if you want epsilon distance for proposed points, set function for epsilon ("line" is default)
  if (opt.state$opt.problem$control$infill.eps.proposed.points.rf) {
    control = opt.state$opt.problem$control
    eps_start = control$infill.eps.start.proposed.points.rf
    if (control$infill.eps.fkt.proposed.points.rf == "line") {
      eps_fkt = function(x) {
        - eps_start / (control$iters - 1) * x + eps_start / (control$iters - 1) * control$iters
      }
    } else if (control$infill.eps.fkt.proposed.points.rf == "parable") {
      # f(x) = a*x^2 + b*x + c
      eps_fkt = function(x) {
        (eps_start / (control$iters^2 - 2 * control$iters + 1)) * x^2 +  
          (-2 * control$iters * eps_start / (control$iters^2 - 2 * control$iters + 1)) * x + 
          (2 * control$iters * eps_start / (control$iters^2 - 2 * control$iters + 1) - 
             eps_start / (control$iters^2 - 2 * control$iters + 1) + eps_start) 
      }
    } else if (control$infill.eps.fkt.proposed.points.rf == "neg.parable") {
      # f(x) = - a*x^2 + b*x + c
      eps_fkt = function(x) {
        - (eps_start - (2 * eps_start * control$iters - eps_start * control$iters^2) / (2 * control$iters - control$iters^2 - 1)) * x^2 + 
          (2 * eps_start - 2 * (2 * eps_start * control$iters - eps_start * control$iters^2) / (2 * control$iters - control$iters^2 -1)) * x + 
          ((2 * eps_start * control$iters - eps_start * control$iters^2) / (2 * control$iters - control$iters^2 -1))
      }
    } else {
      stop("You can only set line, parable or neg.parable for eps.fkt.proposed.points.rf")
    }
    opt.state$opt.problem$control$infill.eps = eps_start
  }
  
  repeat {
    prop = proposePoints.OptState(opt.state)
    
    # if you want epsilon distance for proposed points, save eps in path
    if (opt.state$opt.problem$control$infill.eps.proposed.points.rf) {
      opt.state$opt.path$env$eps[opt.state$loop] = opt.state$opt.problem$control$infill.eps
    }
    
    evalProposedPoints.OptState(opt.state, prop)
    finalizeMboLoop(opt.state)
    terminate = getOptStateTermination(opt.state)
    
    if (terminate$term) {
      break
    }
    
    # if you want epsilon distance for proposed points, lower the epsilon in each iteration
    if (opt.state$opt.problem$control$infill.eps.proposed.points.rf) {
      opt.state$opt.problem$control$infill.eps = eps_fkt(opt.state$loop)
    }
  }
  opt.state
}

finalizeMboLoop = function(opt.state) {
  # put resampling of surrogate learner and the model itself in the result environment
  opt.result = getOptStateOptResult(opt.state)
  setOptResultResampleResults(opt.result, opt.state)
  setOptResultStoredModels(opt.result, opt.state)
  # Indicate, that we are finished by increasing the loop by one
  setOptStateLoop(opt.state)
  # save on disk routine
  # save with increased loop so we can directly start from here again
  if (getOptStateShouldSave(opt.state))
    saveOptState(opt.state)
  invisible()
}
