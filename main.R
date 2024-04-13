library(smoof)
library(vioplot)

prsAlgorithm <- function(func, dimentions, nPoints){
  f <- func(dimentions)
  minVal <- NULL
  
  for(x in 1:nPoints){
    point <- runif(dimentions, getLowerBoxConstraints(f), getUpperBoxConstraints(f))
    value <- f(point)
    
    if(is.null(minVal) || value <  minVal){
      minVal <- value
    }
  }
  
  return(minVal)
}

msAlgorithm <- function(func, dimentions, nPoints){
  f <- func(dimentions)
  counter <- 0
  minVal <- NULL
  
  for(x in 1:nPoints){
    point <- runif(dimentions, getLowerBoxConstraints(f), getUpperBoxConstraints(f))
    
    result <- optim(point, f, method = "L-BFGS-B", lower = getLowerBoxConstraints(f), upper = getUpperBoxConstraints(f))
    value <- as.numeric(result$value)
    counter <- counter + as.numeric(result$counts[1])
    
    if(is.null(minVal) || value < minVal){
      minVal <- value
    }
  }
  return(list(minVal, counter))
}

compareAlgorithms <- function(alg1, alg2, func, dimentions){
  alg1Result <- replicate(50, alg1(func, dimentions, 100))
  alg1Average <- as.numeric(alg1Result[2,])
  alg1Points <- as.numeric(alg1Result[1,])
  average <- mean(alg1Average)
  alg2Result <- replicate(50, alg2(func, dimentions, average))
  
  return(list(alg1Points, alg2Result))
}

rastrigin2Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRastriginFunction, 2)
hist(rastrigin2Result[[1]], main = "Histogram algorytmu MS (Rastrigin, 2D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rastrigin2Result[[2]], main = "Histogram algorytmu PRS (Rastrigin, 2D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rastrigin2Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rastrigin, 2D)", ylab = "Zakres minimów")
vioplot(rastrigin2Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rastrigin, 2D)", ylab = "Zakres minimów")
t.test(x = rastrigin2Result[[1]], conf.level = 0.95)
t.test(x = rastrigin2Result[[2]], conf.level = 0.95)
t.test(x = rastrigin2Result[[2]], y = rastrigin2Result[[1]], conf.level = 0.95)

rosenbrock2Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRosenbrockFunction, 2)
hist(rosenbrock2Result[[1]], main = "Histogram algorytmu MS (Rosenbrock, 2D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rosenbrock2Result[[2]], main = "Histogram algorytmu PRS (Rosenbrock, 2D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rosenbrock2Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rosenbrock, 2D)", ylab = "Zakres minimów")
vioplot(rosenbrock2Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rosenbrock, 2D)", ylab = "Zakres minimów")
t.test(x = rosenbrock2Result[[1]], conf.level = 0.95)
t.test(x = rosenbrock2Result[[2]], conf.level = 0.95)
t.test(x = rosenbrock2Result[[2]], y = rosenbrock2Result[[1]], conf.level = 0.95)

rastrigin10Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRastriginFunction, 10)
hist(rastrigin10Result[[1]], main = "Histogram algorytmu MS (Rastrigin, 10D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rastrigin10Result[[2]], main = "Histogram algorytmu PRS (Rastrigin, 10D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rastrigin10Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rastrigin, 10D)", ylab = "Zakres minimów")
vioplot(rastrigin10Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rastrigin, 10D)", ylab = "Zakres minimów")
t.test(x = rastrigin10Result[[1]], conf.level = 0.95)
t.test(x = rastrigin10Result[[2]], conf.level = 0.95)
t.test(x = rastrigin10Result[[2]], y = rastrigin10Result[[1]], conf.level = 0.95)

rosenbrock10Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRosenbrockFunction, 10)
hist(rosenbrock10Result[[1]], main = "Histogram algorytmu MS (Rosenbrock, 10D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rosenbrock10Result[[2]], main = "Histogram algorytmu PRS (Rosenbrock, 10D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rosenbrock10Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rosenbrock, 10D)", ylab = "Zakres minimów")
vioplot(rosenbrock10Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rosenbrock, 10D)", ylab = "Zakres minimów")
t.test(x = rosenbrock10Result[[1]], conf.level = 0.95)
t.test(x = rosenbrock10Result[[2]], conf.level = 0.95)
t.test(x = rosenbrock10Result[[2]], y = rosenbrock10Result[[1]], conf.level = 0.95)

rastrigin20Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRastriginFunction, 20)
hist(rastrigin20Result[[1]], main = "Histogram algorytmu MS (Rastrigin, 20D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rastrigin20Result[[2]], main = "Histogram algorytmu PRS (Rastrigin, 20D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rastrigin20Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rastrigin, 20D)", ylab = "Zakres minimów")
vioplot(rastrigin20Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rastrigin, 20D)", ylab = "Zakres minimów")
t.test(x = rastrigin20Result[[1]], conf.level = 0.95)
t.test(x = rastrigin20Result[[2]], conf.level = 0.95)
t.test(x = rastrigin20Result[[2]], y = rastrigin20Result[[1]], conf.level = 0.95)

rosenbrock20Result <- compareAlgorithms(msAlgorithm, prsAlgorithm, makeRosenbrockFunction, 20)
hist(rosenbrock20Result[[1]], main = "Histogram algorytmu MS (Rosenbrock, 20D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
hist(rosenbrock20Result[[2]], main = "Histogram algorytmu PRS (Rosenbrock, 20D)", xlab = "Zakresy minimów", ylab = "Ilość minimów")
vioplot(rosenbrock20Result[[1]], main = "Wykres skrzypcowy algorytmu MS (Rosenbrock, 20D)", ylab = "Zakres minimów")
vioplot(rosenbrock20Result[[2]], main = "Wykres skrzypcowy algorytmu PRS (Rosenbrock, 20D)", ylab = "Zakres minimów")
t.test(x = rosenbrock20Result[[1]], conf.level = 0.95)
t.test(x = rosenbrock20Result[[2]], conf.level = 0.95)
t.test(x = rosenbrock20Result[[2]], y = rosenbrock20Result[[1]], conf.level = 0.95)
