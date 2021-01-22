# Generates Erdos-Renyi Random Graph
ER <- function(p, beta){
   # p choose 2 possible edges e without self-loops
   e <- (p*(p-1))/2
   # make adjacency matrix X
   X <- matrix(0, p, p)
   # create e possible random edges, X is a Bernoulli(beta) random variable
   edges <- rbinom(e, 1, beta)
   # Fill X
   X[upper.tri(X, diag = FALSE)] <- edges
   # make X symmetric
   X[lower.tri(X)] = t(X)[lower.tri(X)]
   return(X)
}

# Generates the connection probability block
conn <- function(p, beta){
   # p^2 possible ones and zeros
   e <- p^2
   # create e possible random edges, X is a Bernoulli(beta) random variable
   edges <- rbinom(e, 1, beta)
   # Make matrix
   X <- matrix(edges, nrow = p, ncol = p)
   return(X)
}

# Generates Stochastic Block Model matrix
SBM <- function(p, beta_11, beta_12, beta_22){
   if(p%%2 == 0){
   v <- p/2
   
   A_11 <- ER(v, beta_11)
   A_12 <- conn(v, beta_12)
   top <- cbind(A_11, A_12)
   A_22 <- ER(v, beta_22)
   A_21 <- t(A_12)
   bottom <- cbind(A_21, A_22)
   A <- rbind(top, bottom)
   A
   }
}

# Minimizes the quadratic form of the Laplacian matrix
conductance <- function(A){
   p <- ncol(A)
   # find vertex degrees
   d <- colSums(A)
   # find graph laplacian matrix
   L <- diag(d) - A
   # find eigendecomposition of L
   e <- eigen(L)
   # We minimize e subject to the eigenvector. 
   # The smallest eigenvector had eigenvalue 0 so
   # pick the second smallest
   v <- e$vectors[,p-1]
   # component subset assignment
   1*(v>0)
   # Recall that the corresponding eigenvalue is a 
   # metric of the quality of the cut
}

# Maximize the quadratic form of the modularity matrix
mod <- function(A){
   # find vertex degrees
   d <- colSums(A)
   # sum of the degrees
   # (expected number of edges)
   twon <- sum(d)
   # B = v^{T} * (1/2m * (adjacency matrix - dd^{T}/2m)) * v
   B <- A - outer(d, d, "*")/twon
   e <- eigen(B)
   # we pick the dominant eigenvector to maximize
   v <- e$vectors[,1]
   # component subset assignment
   1*(v>0)
}

# clique percolation method
# where k is the minimum size of a clique

# Steps:
# S1: All k-cliques present in A are extracted.
# S2: A new graph, the clique-graph, cg is formed
# where each node represents an identified clique.
# and 2 vertices in cg are connected by an edge only
# if they have k-1 common vertices. Create graph object CG.
# S3: Connected components in the CG are identified.
# S4: Each connected component in CG represents a 
# community.
# S5: comm is the set of communities formed for CG. 
CPM <- function(A, k=3){
   g <- graph.adjacency(A, mode = "undirected")
   cliq <- max_cliques(g, min=k)
   # number of cliques
   n_cliq <- length(cliq)
   # clique graph
   cg <- matrix(0 , n_cliq, n_cliq)
   for(i in 1:(n_cliq - 1)) {
      for(j in (i+1):n_cliq){
         if(length(intersect(cliq[[i]], cliq[[j]])) >= k-1){
            cg[i,j] <- 1
            cg[j,i] <- 1
         }
      }
   }
   
   # clique graph obj
   CG <- graph.adjacency(cg, mode = "undirected")
   # connected components
   comp <- components(CG)$membership
   # communities
   comm <- vector('list', max(comp))
   for( j in 1:(max(comp))){
      for( i in which(comp == j)){
         comm[[j]] <- union(comm[[j]], cliq[[i]])
      }
   }
   comm
}  

# This vertex coloring works when the communities are linearly separable
# but if not then the color assignments get overwritten.

#   sep <- rep(0, ncol(A))
#   ct <- 1
#   for(vec in comm){
#      for(i in 1:length(vec)){
#         sep[vec[i]] <- ct
#      }
#      ct <- ct + 1
#   }
#   return(list("comm" = comm, "sep" = sep))
#}

# Degree-corrected Stochastic Block Model Optimization
quality <- function(kap, ems, nc=length(kap)){
   tot <- 0
   for(i in 1:nc){
      for(j in 1:nc){
         if(ems[i, j]>0){
            tot <- tot+ems[i, j]*log(ems[i, j]/(kap[i]*kap[j]))
         }
      }
   }
   tot
}

kapems <- function(b, s, nc=2){
   d <- apply(b, 2, sum)
   kap <- rep(0,nc)
   ems <- matrix(0, nc, nc)
   for(i in 1:nc){
      ip <- which(s==i)
      kap[i] <- sum(d[ip])
      for(j in 1:nc){
         jp <- which(s==j)
         if(length(ip)==1 || length(jp)==1){
            ems[i, j] <- sum(b[ip, jp])
            }else{
               ems[i, j] <- sum(apply(b[ip,jp], 2, sum))
            }
      }
   }
   list(kap=kap, ems=ems)
}

opt.sbm <- function(b, nc=2, T=1e4, gam=.985){
   p <- dim(b)[1]
   t <- 0
   s <- sample.int(nc, p, replace=TRUE)
   qold <- -Inf
   while(t<T){
      pos <- sample.int(p, 1)
      lab <- sample(s[which(b[pos, ]>0)], 1)
      sto <- s[pos]
      s[pos] <- lab
      km <- kapems(b, s, nc)
      qnew <- quality(km$kap, km$ems, nc)
      print(c(qold, qnew))
      if(runif(1)<exp(-(qold-qnew)/(T*(gam^t)))){
         qold <- qnew
         }else{
            s[pos] <- sto
         }
      t <- t+1
      print(t)
   }
   s
}

library(igraph)

A <- SBM(p = 10,
         beta_11 = 0.6,
         beta_22 = 0.8,
         beta_12 = 0.4)
g <- graph.adjacency(A, mode = "undirected")
xy <- tk_coords(tkplot(g))

s <- conductance(A)
plot(g, vertex.color = ifelse(s, "red", "lightblue"), edge.width = 4)

q <- mod(A)
plot(g, vertex.color = ifelse(q, "red", "lightblue"), edge.width = 4)


B <- SBM(p = 80,
         beta_11 = 0.6,
         beta_22 = 0.3,
         beta_12 = 0.4)
c = 5
cpm.B <- CPM(B, k = c)

num.comm <- length(cpm.B)
cmap <- rainbow(num.comm)

f <- graph.adjacency(B, mode = "undirected")
plot(f, vertex.color = "white", edge.width = 1, mark.groups = cpm.B)
title(main = paste0("k: ", c, ", p: ", 80))

opt <- opt.sbm(B, nc=2, T=1e4, gam=.985)
plot(f, vertex.color = cmap[opt], edge.width = 2)


###############################################################################

# round(runif(1),2)
rand1 <- 0.7
rand2 <- 0.7
rand3 <- 0.1

B <- SBM(p = 80,
         beta_11 = rand1,
         beta_22 = rand2,
         beta_12 = rand3)
c = 12
cpm.B <- CPM(B, k = c)

num.comm <- length(cpm.B$comm)
cmap <- rainbow(num.comm)

f <- graph.adjacency(B, mode = "undirected")
plot(f, vertex.color = cmap[cpm.B$sep], edge.width = 1)
title(main = paste0("Beta_11: ", rand1, ", Beta_22: ", rand2, ", Beta_12: ", 
                    rand3, ", k: ", c, ", p: ", 80))




# Generates ordered community
p <- 150
n <- p/5

A_11 <- ER(n, 0.9)
A_12 <- conn(n, 0.2)
A_13 <- matrix(0, n, n)
A_13 <- matrix(0, n, n)
A_14 <- matrix(0, n, n)
A_15 <- matrix(0, n, n)
top <- cbind(A_11, A_12, A_13, A_14, A_15)
A_21 <- matrix(0, n, n)
A_22 <- ER(n, 0.9)
A_23 <- conn(n, 0.2)
A_24 <- matrix(0, n, n)
A_25 <- matrix(0, n, n)
mid2 <- cbind(A_21, A_22, A_23, A_24, A_25)
A_31 <- matrix(0, n, n)
A_32 <- matrix(0, n, n)
A_33 <- ER(n, 0.9)
A_34 <- conn(n, 0.2)
A_35 <- matrix(0, n, n)
mid3 <- cbind(A_31, A_32, A_33, A_34, A_35)
A_41 <- matrix(0, n, n)
A_42 <- matrix(0, n, n)
A_43 <- matrix(0, n, n)
A_44 <- ER(n, 0.9)
A_45 <- conn(n, 0.2)
mid4 <- cbind(A_41, A_42, A_43, A_44, A_45)
A_51 <- matrix(0, n, n)
A_52 <- matrix(0, n, n)
A_53 <- matrix(0, n, n)
A_54 <- matrix(0, n, n)
A_55 <- ER(n, 0.9)
bottom <- cbind(A_51, A_52, A_53, A_54, A_55)
A <- rbind(top, mid2, mid3, mid4, bottom)

cpm.A <- CPM(A, k = 8)

num.comm <- length(cpm.A)
cmap <- rainbow(num.comm)

h <- graph.adjacency(A, mode = "undirected")
plot(h, mark.groups = cpm.A, vertex.color = "white", edge.width = 1)


# Generates assortative community
p <- 150
n <- p/5

A_11 <- ER(n, 0.9)
A_12 <- conn(n, 0.1)
A_13 <- conn(n, 0.1)
A_13 <- conn(n, 0.1)
A_14 <- conn(n, 0.1)
A_15 <- conn(n, 0.1)
top <- cbind(A_11, A_12, A_13, A_14, A_15)
A_21 <- conn(n, 0.1)
A_22 <- ER(n, 0.9)
A_23 <- conn(n, 0.1)
A_24 <- conn(n, 0.1)
A_25 <- conn(n, 0.1)
mid2 <- cbind(A_21, A_22, A_23, A_24, A_25)
A_31 <- conn(n, 0.1)
A_32 <- conn(n, 0.1)
A_33 <- ER(n, 0.9)
A_34 <- conn(n, 0.1)
A_35 <- conn(n, 0.1)
mid3 <- cbind(A_31, A_32, A_33, A_34, A_35)
A_41 <- conn(n, 0.1)
A_42 <- conn(n, 0.1)
A_43 <- conn(n, 0.1)
A_44 <- ER(n, 0.9)
A_45 <- conn(n, 0.1)
mid4 <- cbind(A_41, A_42, A_43, A_44, A_45)
A_51 <- conn(n, 0.1)
A_52 <- conn(n, 0.1)
A_53 <- conn(n, 0.1)
A_54 <- conn(n, 0.1)
A_55 <- ER(n, 0.9)
bottom <- cbind(A_51, A_52, A_53, A_54, A_55)
A <- rbind(top, mid2, mid3, mid4, bottom)

cpm.A <- CPM(A, k = 8)

num.comm <- length(cpm.A)
cmap <- rainbow(num.comm)

h <- graph.adjacency(A, mode = "undirected")
plot(h, mark.groups = cpm.A, vertex.color = "white", edge.width = 1)
