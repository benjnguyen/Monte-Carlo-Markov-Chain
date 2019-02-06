library(lattice)

# Set seed to have reproducible results
set.seed(42)

# Generate configuration (as vector and matrix)
getConfig <- function(nrow, ncol, statespace = c(-1, 1))
{
  config <- sample(x = statespace, size = nrow*ncol, replace = TRUE)
  mat <- matrix(config, nrow = nrow, byrow = TRUE)
  result <- list(vec = config, mat = mat)
  return(result)
}

# Find right neighbor of a selected vertex
right <- function(x)
{
  r <- x
  len <- length(x)
  for (i in 1:len)
  {
    if (i != len)r[i] <- x[i+1]
    if (i == len)r[i] <- x[1]
  }
  return(r)
}

# Find left neighbor of a selected vertex
left <- function(x)
{
  l <- x
  len <- length(x)
  for (i in 1:len)
  {
    if (i != 1)l[i] <- x[i-1]
    if (i == 1)l[i] <- x[len]
  }
  return(l)
}

# Main function simulation
wolffsim <- function(n1, n2, temp, nsteps)
{
  
  ## Initializing parameters and neighbors
  n1 = n1 #nrow
  n2 = n2 #ncol
  N = n1 * n2 #size of lattice
  Temp = temp
  B = 1/Temp
  J = 1
  probAdd = 1 - exp(-2*B*J)
  nsteps = nsteps
  # Initialize configuration to "S" for starting configuration
  S <- getConfig(n1, n2)
  Svec <- S$vec
  Smat <- S$mat
  IndexMat <- matrix(1:N, nrow = n1, ncol = n2, byrow = TRUE)
  # Compute neighbors of every point on the lattice
  nbr = list()
  for (i in 1:N)
  {
    nbr[[i]] = c(u = (i - n2 - 1) %% N + 1, 
                 d = (i + n2 - 1) %% N + 1,
                 r = 0,
                 l = 0)
  }
  l = apply(IndexMat, 1, left)
  r = apply(IndexMat, 1, right)
  for(i in 1:length(nbr))
  {
    nbr[[i]][3] <- r[i]
    nbr[[i]][4] <- l[i]
  }
  
  
  ## Initializing repository for data
  SaveMC <- list()
  SaveMC[[1]] <- Smat
  timestep <- c(1:nsteps)
  magnetization <- c(1:nsteps)
  
  ## Running the algorithm
  for (step in 1:nsteps)
  {
    k = sample(1:N, 1)
    Pocket = list(k)
    Cluster = list(k)
    while (length(Pocket) != 0)
    {
      j = sample(Pocket, 1)
      for (l in nbr[[unlist(j)]])
      {
        if (Svec[l] == Svec[unlist(j)])
        {
          if (runif(1) < probAdd)
          {
            if (!(l %in% Cluster))
            {
              Pocket <- append(Pocket, list(l))
              Cluster <- append(Cluster, list(l))
            }
          }
        }
      }
      Pocket <- Pocket[Pocket != unlist(j)]
      #print(length(Cluster))
    }
    for (t in Cluster)
    {
      Svec[t] <- -1 * Svec[t]
    }
    print(step)
    timestep[step] <- step
    magnetization[step] <- mean(Svec)
    SaveMC[[step+1]] <- matrix(Svec, nrow = n1, ncol = n2)
  }
  return(list(timestep = timestep, mag = magnetization, MC = SaveMC))
}

# Various cases

# 50x50 lattice, temp 0.5, 100 iterations
test <- wolffsim(50, 50, 0.5, 100)
ttime <- test$timestep
tmag <- test$mag

# 50x50 lattice, temp 2.26, 100 iterations
test1 <- wolffsim(50, 50, 2.26, 100)
t1time <- test1$timestep
t1mag <- test1$mag

# 50x50 lattice, temp 10, 100 iterations
test2 <- wolffsim(50, 50, 10, 100)
t2time <- test2$timestep
t2mag <- test2$mag

# Code example for generating plots with saved time and magnetization values
plot(ttime, tmag,
     col = "red",
     xlab = "Wolff Sweeps",
     ylab = "Total Magnetization",
     main = "Using Total Magnetization to Observe Critical Slowing Down")
lines(ttime, tmag, col = "red", lwd = 0.5, lty = 1)
points(t1time, t1mag, col = "blue")
lines(t1time, t1mag, col= "blue", lwd = 2)
points(t2time, t2mag, col = "green")
lines(t2time, t2mag, col = "green", lwd = 2)
legend(3, 0.75, legend=c("Temp = 0.5","Temp = 2.26", "Temp = 10"),
       col=c("red", "blue", "green"), pch = 21)


# Sometimes, we want to continue the algorithm given any lattice
# This function does not generate a new configuration; it uses a user inputted configuration
# and continues the simulation
continuewolffsim <- function(mat, n1, n2, temp, nsteps)
{
  ## Initializing parameters and neighbors
  n1 = n1 #nrow
  n2 = n2 #ncol
  N = n1 * n2
  Temp = temp
  B = 1/Temp
  J = 1
  probAdd = 1 - exp(-2*B*J)
  nsteps = nsteps
  #S <- getConfig(n1, n2)
  S <- mat
  Svec <- S$vec
  Smat <- S$mat
  IndexMat <- matrix(1:N, nrow = n1, ncol = n2, byrow = TRUE)
  nbr = list()
  for (i in 1:N)
  {
    nbr[[i]] = c(u = (i - n2 - 1) %% N + 1, 
                 d = (i + n2 - 1) %% N + 1,
                 r = 0,
                 l = 0)
  }
  l = apply(IndexMat, 1, left)
  r = apply(IndexMat, 1, right)
  for(i in 1:length(nbr))
  {
    nbr[[i]][3] <- r[i]
    nbr[[i]][4] <- l[i]
  }
  
  
  ## Initializing repository for data
  SaveMC <- list()
  SaveMC[[1]] <- Smat
  timestep <- c(1:nsteps)
  magnetization <- c(1:nsteps)
  
  ## Running the algorithm
  for (step in 1:nsteps)
  {
    k = sample(1:N, 1)
    Pocket = list(k)
    Cluster = list(k)
    while (length(Pocket) != 0)
    {
      j = sample(Pocket, 1)
      for (l in nbr[[unlist(j)]])
      {
        if (Svec[l] == Svec[unlist(j)])
        {
          if (runif(1) < probAdd)
          {
            if (!(l %in% Cluster))
            {
              Pocket <- append(Pocket, list(l))
              Cluster <- append(Cluster, list(l))
            }
          }
        }
      }
      Pocket <- Pocket[Pocket != unlist(j)]
      #print(length(Cluster))
    }
    for (t in Cluster)
    {
      Svec[t] <- -1 * Svec[t]
    }
    print(step)
    timestep[step] <- step
    magnetization[step] <- mean(Svec)
    SaveMC[[step+1]] <- matrix(Svec, nrow = n1, ncol = n2)
  }
  return(list(timestep = timestep, mag = magnetization, MC = SaveMC))
}
