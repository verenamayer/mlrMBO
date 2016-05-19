# random Search

randomSearch = function(fun, control, show.info = getOption("mlrMBO.show.info", TRUE), more.args = list()) {
  design = generateRandomDesign(n = 1, par.set = smoof::getParamSet(fun))
  control$infill.crit = "random"
  control$multipoint.cl.lie = min
  mbo(fun = fun, design = design, learner = NULL, control = control, show.info = show.info, more.args = more.args)
}