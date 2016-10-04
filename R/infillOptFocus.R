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

    # do iterations where we focus the region-of-interest around the current best point
    for (local.iter in seq_len(control$infill.opt.focussearch.maxit)) {
      # predict on design where NAs were imputed, but return proposed points with NAs
      newdesign = generateDesign(control$infill.opt.focussearch.points, ps.local, randomLHS)
      
      
      # if you want epsilon distance for proposed points, delete points in the newdesign with distance 
      # smaller than epsilon; if all points in the newdesign smaller take a random point.
      if (control$infill.eps.proposed.points.rf) {
        # for gower dist, the factor levels of design and newdesign have to be equal
        for (i in 1:length(newdesign)) {
          if (is.factor(newdesign[, i])) {
            colname = colnames(newdesign)[i]
            if (!length(levels(design[ , colname])) == length(levels(newdesign[ , colname]))) {
              levels(newdesign[ , colname]) = c(levels(newdesign[ , colname]), 
                                                levels(design[ , colname])[which(!(levels(design[ , colname]) %in% 
                                                                                     levels(newdesign[ , colname])))])
              #right order
              newdesign[ , colname] = factor(newdesign[, colname], levels = levels(design[ , colname]))
            }
          }
        }
        # compute distance between designpoints (newdesign to design --> in each row are the distances to a single new point)
        if (control$infill.eps.dist.proposed.points == "gower") {
          ranges = c()
          for (l in seq(1:length(newdesign))) {
            if (is.numeric(design[, l])) {
              minimum = min(design[, l], newdesign[, l])
              maximum = max(design[, l], newdesign[, l])
              ranges[l] = abs(maximum - minimum)
            } else {ranges[l] = 1}  
          }
          dist_designs = proxy::dist(newdesign, design[ , which(colnames(design) != "y"), drop = FALSE], by_rows = TRUE,
                                     method = "gower", ranges = ranges)
          if (is.null(ranges)) {
            warning("ranges = NULL, the default range calculation in proxy::dist isn't like at the original Gower formula! 
                    (it should be over both designs)")
          }
        } else if (control$infill.eps.dist.proposed.points == "euclidean") {
          dist_designs = proxy::dist(newdesign, design[ , which(colnames(design) != "y"), drop = FALSE], by_rows = TRUE,
                                     method = "euclidean")
        }
        
        # if distance < epsilon, delete point in the newdesign
        # in each row of dist_designs: 
        # if at least one value is smaller than eps, the distance of the new point 
        # is too small to at least one design point --> delete new point
        newpoints_dist_smaller_eps = unique(which(dist_designs < control$infill.eps, arr.ind = TRUE)[, 1])
        if (length(newpoints_dist_smaller_eps) != 0 &&
            length(newpoints_dist_smaller_eps) != dim(newdesign)[1]) {
          newdesign = newdesign[-newpoints_dist_smaller_eps, , drop = FALSE]
        
        # if all points < epsilon, take the point with max distance
        # } else if (length(newpoints_dist_smaller_eps) == dim(newdesign)[1]) {
        #   newpoints_dist_max = which(dist_designs == max(dist_designs), arr.ind = TRUE)[, 1]
          # sample one max point; but maybe its better to take the point with the best infill.crit value 
          #if (length(newpoints_dist_max) > 1) {
          #  newpoints_dist_max = sample(newpoints_dist_max, 1)
          #}
        # newdesign = newdesign[newpoints_dist_max, , drop = FALSE]
        # }
        
        # if all points < epsilon, take a random point
        } else if (length(newpoints_dist_smaller_eps) == dim(newdesign)[1]) {
          random.point = sample(newpoints_dist_smaller_eps, 1)
          newdesign = newdesign[random.point, , drop = FALSE]
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
               par$values = par$values[-to.del, drop = FALSE]
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
