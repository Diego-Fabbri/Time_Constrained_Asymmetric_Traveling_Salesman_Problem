#Set your own working directory
setwd("C:/Users/diego/Documents/R/Projects/GitHub_Projects/Optimization/Asymmetric Traveling Salesman Problem with Time Windows")

# Import lpSolve package
library(lpSolve)

#Import required packages (ompr)
library(dplyr)
library(ROI)
library(ROI.plugin.symphony)
library(ompr)
library(ompr.roi)

#Set t0
t0 <- 0

#Set speed
speed <- 1

#Set matrix of costs (distances)
c <- matrix(c(10000, 5.5, 4.2, 2.6, 2.4, 1.3, 2.5, 4.3,
              4.7, 10000, 3.7, 2.1, 5.1, 6, 7.2, 9,
              4.2, 4.5, 10000, 1.6, 3.2, 5.5, 6.7, 8.5,
              2.6, 2.9, 1.6, 10000, 3, 3.9, 5.1, 6.9,
              3.8, 4.1, 2.8, 1.2, 10000, 5.1, 6.3, 8.1,
              3.9, 7.4, 6.1, 4.5, 3.3, 10000, 1.2, 3,
              3.5, 7, 5.7, 4.1, 2.9, 1.2, 10000, 2.3,
              5.8, 9.3, 8, 6.4, 5.2, 3, 2.3, 10000), nrow = 8, byrow = TRUE)



#Set n (depot + nodes to visit = problem size)
n <- ncol(c)

#Set matrix of travel duration
travel_time <- array(dim = c(n, n))

for (i in 1:n) {
  for (j in 1:n) {
    travel_time[i, j] <- c[i, j]/speed
  }
}

#Set Lower Bounds
LB <- c(0, 0, 0, 0, 0, 0, 0, 0)

#Set Upper Bounds
UB <- c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000)

#Set matrix M
M <- array(dim = c(n, n))

for (i in 1:n) {
  for (j in 1:n) {
    M[i, j] <- UB[i] - LB[i] + travel_time[i, j]
  }
}

#Build Model
Model <- MIPModel() %>%
  add_variable(x[i], i = 1:(n+1), type = "continuous") %>% #define variables
  add_variable(y[i,j], i = 1:(n+1), j = 1:(n+1), type = "binary") %>%
  set_objective(x[n+1], "min") %>% #define objective
  add_constraint(x[1] == t0) %>% #define constraints
  add_constraint(y[i, i] == 0, i = 1:(n+1)) %>%
  add_constraint(x[i] >= travel_time[1,i]*y[1,i], i = 2:n) %>% 
  add_constraint(x[i] - x[j] <= UB[i] - LB[j] -M[i, j]*y[i, j], i = 2:n, j = 2:n, i!=j) %>%
  add_constraint(sum_expr(y[i, j], i = 1:n) == 1, j = 2:n) %>%
  add_constraint(sum_expr(y[i, j], j = 1:n) == 1, i = 2:n) %>%
  add_constraint(x[i] + travel_time[i,1] <= x[n+1], i = 2:n) %>%
  add_constraint(x[i] >= LB[i], i = 2:n) %>%
  add_constraint(x[i] <= UB[i], i = 2:n) %>%
  add_constraint(sum_expr(y[i, 1], i = 2:n) == 1) %>%
  add_constraint(sum_expr(y[1, j], j = 2:n) == 1) %>%
  solve_model(with_ROI(solver = "symphony", verbosity = 1))

#Model summary
##Status
print(paste("Model status is:", Model$status))

##Objective Function
print(paste("Objective value:", objective_value(Model)))


#x variables 
for (a in 1:(n+1)) {
  tmp_x <- get_solution(Model, x[i]) %>%
    filter(variable == "x", i == a) %>%
    select(value)
  
  if (tmp_x != 0) {
    print(paste("--->x[", (a-1), "] =", tmp_x))
  }
}

#y variables
for (a in 1:(n+1)) {
  for (b in 1:(n+1)) {
    if(a!=b){
      tmp <- get_solution(Model, y[i,j]) %>%
        filter(variable == "y", i == a , j == b) %>%
        select(value)
      if (tmp != 0) {
        print(paste("--->y[", a,",",b, "] =", tmp))
      }
    }
  }
}


