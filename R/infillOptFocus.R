# Random search, where we shrink the region of interest after restarts
# around the currently best point. only numeric / ints are currently "shrunken"
# works for ALL parameter sets
#
#FIXME it would be nice to have a REASONABLE way to shrink categorical stuff too.
#FIXME should we shrink if a local value is NA (dependent param)
#
# See infillOptCMAES.R for interface explanation.
infillOptFocus = function(infill.crit, models, control, par.set, opt.path, design, iter, ...) {
  global.y = Inf
  
  # restart the whole crap some times
  for (restart.iter in seq_len(control$infill.opt.restarts)) {
    # copy parset so we can shrink it
    ps.local = par.set
    if (control$infill.eps.proposed.points.rf == TRUE) {
      eps = control$infill.eps
    }

    # do iterations where we focus the region-of-interest around the current best point
    for (local.iter in seq_len(control$infill.opt.focussearch.maxit)) {
      # predict on design where NAs were imputed, but return proposed points with NAs
      newdesign = generateDesign(control$infill.opt.focussearch.points, ps.local, randomLHS)
      
      # if you want epsilon distance for proposed points, delete points in the newdesign with distance 
      # smaller than epsilon; if all points in the newdesign smaller take the point with max distance.
      if (control$infill.eps.proposed.points.rf == TRUE) {
        # compute distance between designpoints and save points with dist < epsilon
        dist_designs = gower.dist(data.x = newdesign, data.y = design[,-dim(design)[2], drop = FALSE])
        newpoints_dist_smaller_eps = unique(which(dist_designs < eps, arr.ind = TRUE)[, 1])
        
        # if distance < epsilon, delete point in the newdesign
        if (length(newpoints_dist_smaller_eps) != 0 &&
            length(newpoints_dist_smaller_eps) != dim(newdesign)[1]) {
          newdesign = newdesign[-newpoints_dist_smaller_eps, , drop = FALSE]
        }
        
        # if all points < epsilon, take the point with max distance
        if (length(newpoints_dist_smaller_eps) == dim(newdesign)[1]) {
          newpoints_dist_max = which(dist_designs == max(dist_designs), arr.ind = TRUE)[, 1]
          # sample one point; but maybe its better to take the point with the best infill.crit value 
          #if (length(row_dist_max) > 1) {
          #  row_dist_max = sample(row_dist_max, 1)
          #}
          newdesign = newdesign[newpoints_dist_max, , drop = FALSE]
        }
      }
      
      # convert to param encoding our model was trained on and can use
      newdesign = convertDataFrameCols(newdesign, ints.as.num = TRUE, logicals.as.factor = TRUE)
      y = infill.crit(newdesign, models, control, ps.local, design, iter, ...)

      # get current best value
      local.index = getMinIndex(y, ties.method = "random")
      local.y = y[local.index]
      local.x.df = newdesign[local.index, , drop = FALSE]
      local.x.list = dfRowToList(recodeTypes(local.x.df, ps.local), ps.local, 1)

      # if we found a new best value, store it
      if (local.y < global.y) {
        global.x.df = local.x.df
        global.y = local.y
      }

      # now shrink ps.local object so we search more locally
       ps.local$pars = lapply(ps.local$pars, function(par) {
         # only shrink when there is a value
         val = local.x.list[[par$id]]
         if (!isScalarNA(val)) {
           if (isNumeric(par)) {
             # shrink to range / 2, centered at val
             range = par$upper - par$lower
             par$lower = pmax(par$lower, val - (range / 4))
             par$upper = pmin(par$upper, val + (range / 4))
             if (isInteger(par)) {
               par$lower = floor(par$lower)
               par$upper = ceiling(par$upper)
             }
           } else if (isDiscrete(par)) {
             # randomly drop a level, which is not val
             if (length(par$values) > 1L) {
               val.names = names(par$values)
               # remove current val from delete options, should work also for NA
               val.names = setdiff(val.names, val)
               to.del = sample(seq_along(val.names), 1)
               par$values = par$values[-to.del]
             }
           }
         }
         return(par)
      })
    }
  }
  recodeTypes(global.x.df, par.set)
}

# as we operate on other types for the learner (ints are nums, logs are factors),
# we have to convert back sometimes for dfRowsToList to work
recodeTypes = function(df, par.set) {
  types = unlist(lapply(par.set$pars, function(p) rep(p$type, p$len)))
  for (i in seq_col(df)) {
    if (types[i] %in% c("integer", "integervector"))
      df[, i] = as.integer(df[,i])
    else if (types[i] %in% c("logical", "logicalvector"))
      df[, i] = as.logical(as.character(df[,i]))
  }
  return(df)
}
