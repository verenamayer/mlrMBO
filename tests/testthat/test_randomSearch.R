context("randomSearch")

test_that("basic randomSearch works", {
  ctrl = makeMBOControl()
  ctrl = setMBOControlTermination(ctrl, iters = 20L)
  or = randomSearch(fun = testf.fsphere.2d, control = ctrl, show.info = FALSE)
  expect_true(!is.na(or$y))
  expect_equal(or$y, testf.fsphere.2d(or$x))
  expect_equal(getOptPathLength(or$opt.path), 20)
  expect_true(is.list(or$x))
  expect_equal(names(or$x), names(testp.fsphere.2d$pars))
  
  #FIXME, check result if time is met before iters
  objfun = smoof::makeSingleObjectiveFunction(
    fn = function(x) {
      x = (x$x-5)^2
      attr(x,"exec.time") = 10L
      return(x)
    },
    par.set = makeParamSet(
      makeNumericParam("x", lower = 0, upper = 10)
    ),
    has.simple.signature = FALSE
  )
  ctrl2 = setMBOControlTermination(ctrl, exec.time.budget = 100)
  or = randomSearch(fun = objfun, control = ctrl2, show.info = FALSE)
  expect_equal(sum(as.data.frame(or$opt.path)$exec.time), 150)
  expect_equal(getOptPathLength(or$opt.path), 10)

  #now with parallel execution
  ctrl3 = makeMBOControl(schedule.nodes = 2)
  ctrl3 = setMBOControlTermination(ctrl3, time.budget = 5)
  or = randomSearch(fun = objfun, control = ctrl3, show.info = FALSE)
})

