#--------------------------------------------------------
# Keith Gladstone
# ORF 360 Final Project Code
# May 2016
#--------------------------------------------------------
library(sfsmisc)
#--------------------------------------------------------
# Adjustable parameters
#--------------------------------------------------------
k <- 7 # constant initial capacity on every flight-leg
T <- 30 # number of selling periods

#--------------------------------------------------------
# Process Product Incidence Matrix
# Format: 
# Columns are products (itineraries)
# Rows are resources (flight-legs)
# Second-to-last row is product prices
# Last row is product request frequencies
#--------------------------------------------------------
A <- read.csv('csv/product_incidence.csv', 
              stringsAsFactors = FALSE)
rownames(A) <- A[,1]
A <- A[,2:length(A)]
n <- length(A) # number of products (itineraries)
m <- nrow(A) - 2 # number of resources (flight-legs)
p <- as.vector(A[m + 1,], 
               'numeric') # price of product j 
q <- as.vector(A[m + 2,], 
               'numeric') # arrival freq. for product j
A <- A[1:m,] # remove prices, frequencies from A
x.count <- (k+1)^m # (k+1)^m unique state vectors  
u_alt <- digitsBase(1:(2^(n)-1),2) # decision vectors
V <- matrix(0, nrow = (T + 1), ncol = x.count) 

# IMPORTANT NOTE: Interpretation of V Matrix
# -----
# Columns from left to right: 
# (States x, indexed from 1 to (k+1)^m, 
#   where m is number of resources and k
#   is starting capacity for all flight-legs)
# Rows from top to bottom: (Time T + 1 through 1)
#

#--------------------------------------------------------
# Vector mapping functions
#--------------------------------------------------------

# FUNCTION to map a key to an x-vector
stateGivenIndex <- function(key) {
  if (key == 1)
    return(rep(0,m))
  vec <- as.vector(digitsBase(key - 1, k + 1))
  return(c(rep(0,m - length(vec)), vec))
}

# FUNCTION to map an x-vector to a key
indexOfState <- function(vec) {
  result <- 0
  for (place in 1:length(vec)) 
    result <- result + 
      vec[place]*(k+1)^(length(vec) - place)
  return(result + 1)
}

#--------------------------------------------------------
# Bellman building block functions
#--------------------------------------------------------

# FUNCTION to compute opportunity cost of accepting `j'
opportunity_cost <- function(t, x, j) {
  return(V[t+1, indexOfState(x)] - V[(t+1), 
                              indexOfState(x - A[,j])])
}

# FUNCTION to test if x - Aj nonnegative
isFeasible <- function(x, j) {
  return(sum(x - A[,j] < 0) == 0)
}

# Value-iteration
t_0 <- 1 
for (t in T:t_0) {
  print(paste('Time:', t)) # track execution progress
  for (state.index in 1:x.count) {
    x <- stateGivenIndex(state.index) 
    max_V_alt <- 0 
    for (u.index in 1:ncol(u_alt)) {
      u <- u_alt[, u.index] # select decision vector

      # Perform feasibility check
      u_allowed <- TRUE
      for (j in 1:n) 
        if (u[j] == 1 & isFeasible(x,j) == FALSE) {
            u_allowed = FALSE
            break
          } 

      # If this decision is feasible, consider it for V
      if (u_allowed == TRUE) {
        # Summation
        V_alt = 0
        for (j in 1:n)  
          if (u[j] == 1) 
            V_alt <- V_alt + 
                q[j]*u[j]*
                  (p[j] - opportunity_cost(t, x, j))
        if (V_alt > max_V_alt)  
          max_V_alt <- V_alt
      } 
    }

    V[t, state.index] <- V[t+1, state.index] + max_V_alt
  }
}
print('Value Iteration Complete')
write.csv(V, 'csv/value_iteration_output.csv')
print('The value of V(t = 1, x = full capacity) is:')
print(V[1, x.count])