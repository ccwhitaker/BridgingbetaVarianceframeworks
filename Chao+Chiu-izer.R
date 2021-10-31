#chaochiuizeradapted from Tina Harrison"s work


##################################################
library(tidyverse)
data1985 <-
  read.table("percov1985.txt", header = TRUE) %>% select(-1) %>% as.matrix()
data1983 <-
  read.table("percov1983.txt", header = TRUE) %>% select(-1) %>% as.matrix()

chaochiuizer <- function(Z, q, output) {
  

  ##################################################
  # eqn 1, defining a site-species abundance matrix
  N <- ncol(Z) # number of communities
  S <- nrow(Z) # number of species (gamma)
  # drop all rows and columns that sum to 0
  Z <- Z[rowSums(Z) > 0, colSums(Z) > 0]

  ##################################################
  # eqn 2, matrix of mean abundance (the "null")
  M = matrix(rep(rowSums(Z) / N, each = N), ncol = N, byrow = T)

  ##################################################
  #Variance Framework
  #eqn 3, total centered sum of squares
  SumSqrs = sum((Z - M) ^ 2)

  ##################################################
  #Decompostion Framework

  # eqn 5c, multiplicative beta component
  # first, redefine reality to make life a little easier
  Ln = function(x) {
    ifelse(x != 0, log(x), 0)
  }
  Exp = function(x, y) {
    ifelse(x != 0, x ^ y, 0)
  }

  # eqn 5a, gamma diversity. note relativized by total abundance
  D_gamma = ifelse(q != 1,
                   sum(Exp(rowSums(Z)/sum(Z), q))^(1/(1 - q)),
                   exp(-sum((rowSums(Z)/sum(Z))*Ln(rowSums(Z)/sum(Z)))) )

  # eqn 5b, alpha diversity.
  D_alpha = ifelse(q != 1,
                   (1/N)*sum(Exp(Z/sum(Z), q))^(1/(1-q)), 
                   exp(-sum((Z/sum(Z))*Ln(Z/sum(Z))) - log(N)) )

  # eqn 5c, multiplicative beta component
  D_beta = D_gamma/D_alpha

  ##################################################
  # Decompostion normalized overlap measures
  # eqn 6a, N-community Sorensen; effective average proportion of a community's species that are shared across all communities; 1 - CqN is the effective average of unique species in a community
  C_qN = ifelse(q != 1,
                ((1 / D_beta) ^ (q - 1) - (1 / N) ^ (q - 1)) / (1 - (1 / N) ^ (q - 1)),
                1 - (1 / (sum(Z) * Ln(N))) * sum(Z * Ln(Z / rowMeans(Z))))
  Sorensen = 1 - C_qN
  # when Z is within-community rel abundances, q=1 reduces to Horn and q=2 reduces to Morisita-Horn. Is there a general solution for q=1?

  # eqn 6b, N-community Jaccard; effective proportion of species in pooled community that are shared across all communities; 1 - UqN is the effective proportion of unique species in the pooled community
  U_qN = ifelse(q != 1,
                ((1 / D_beta) ^ (1 - q) - (1 / N) ^ (1 - q)) / (1 - (1 / N) ^ (1 - q)),
                1 - (1 / (sum(Z) * Ln(N))) * sum(Z * Ln(Z / rowMeans(Z))))
  Jaccard = 1 - U_qN
  # when Z is within-community rel abundances, q=1 reduces to Horn and q=2 reduces to regional species-overlap from Chiu, Jost & Chao 2014

  #################################################
  #eqn 7a q-th order divergence

  qdiv = (1 / (q - 1)) * (sum((Z ^ q) - (M ^ q)))
  #What q would bridge the variance framework with the decomposition framework?

  #################################################
  #eqn 12a, normalized divergence 

  normdiv1=(sum((Z^q)-(M^q)))/(((N^q)-N)*colSums(M^q)[1])

  # summary stuff
  N; sum(Z); S; q

  # rowSums(Z) # total abundance of ith species in pooled community
  # rowSums(Z)/N # average abundance of ith species per community
  # colSums(Z) # size of the kth community

  # Output function depending on output= "everything", "qdiv", or "normdiv1"

  if(output=="everything") {
    out <-
      c(SumSqrs,
        D_gamma,
        D_alpha,
        D_beta,
        Sorensen,
        Jaccard,
        qdiv,
        normdiv1,
        normdiv2)
    names(out) <-
      c(
        "SumSqrs",
        "D_gamma",
        "D_alpha",
        "D_beta",
        "Sorensen",
        "Jaccard",
        "qdiv",
        "normdiv1",
        'normdiv2'
      )
    out
  }
  else {
    if (output == "normdiv1") {
      out2 <- c(normdiv1)
      out2
    }
    else {
      if (output == "Jaccard") {
        out3 <- c(Jaccard)
        out3
      }
      else {
        if (output == "qdiv") {
          out4 <- c(qdiv)
          out4
        }
        else{
          out5 <- NA
          out5
        }
      }
    }
  }
}

###############################################################################
chaochiuizer(Z=data1985,q=2,output = "normdiv1")

#what do you notice about normdiv1 and the Jaccard index when q equals 1?
# what about sumsqrs and qdiv when q equals 2?


###############################################################################
#how does this vary with q?
chaochiuizer(Z=data1985,q=0,output = "qdiv")
chaochiuizer(Z=data1985,q=0.5,output = "qdiv")
chaochiuizer(Z=data1985,q=1,output = "qdiv")
chaochiuizer(Z=data1985,q=1.5,output = "qdiv")
chaochiuizer(Z=data1985,q=2,output = "qdiv")

#So this is what I am trying to do with the for loop. Unfortunately it looks like my 
#if else statement messed up the for loop
###############################################################################

#how might we use the jaccard or normdiv measure to compare communities? 

qute <- seq(0, 2, by=0.1)#create a vector of possible q values


normvector1985 <- NULL#create an empty vector
normvector1983 <- NULL

for (i in seq(qute)) {

  normvector1985[i] <-
    chaochiuizer(Z = data1985, output = "normdiv1", q = qute[i])

  normvector1983[i] <-
    chaochiuizer(Z = data1983, output = "normdiv1", q = qute[i])

}

qvector1983 <- NULL
qvector1985 <- NULL

for (i in seq(qute)) {

  qvector1985[i] <-
    chaochiuizer(Z = data1985, output = "qdiv", q = qute[i])

  qvector1983[i] <-
    chaochiuizer(Z = data1983, output = "qdiv", q = qute[i])

}

jaccardvect1983<-NULL
jaccardvect1985<-NULL

for (i in seq(qute)) {
  jaccardvect1983[i] <-
    chaochiuizer(Z = data1985, output = "Jaccard", q = qute[i])
  jaccardvect1985[i] <-
    chaochiuizer(Z = data1983, output = "Jaccard", q = qute[i])
}

#figure 2a
notnormalized_divergence <- data.frame(qute, qvector1983, qvector1985) %>%
  pivot_longer(2:3, names_to = "year", values_to = "differentation") %>%
  ggplot(aes(
    x = qute,
    y = differentation,
    group = year,
    color = year
  )) + geom_line() + geom_point()
notnormalized_divergence

#figure 2b
normalized_divergence <-
  data.frame(qute, normvector1983, normvector1985) %>%
  pivot_longer(2:3, names_to = "year", values_to = "differentation") %>%
  ggplot(aes(
    x = qute,
    y = differentation,
    group = year,
    color = year
  )) + geom_line() + geom_point()
normalized_divergence

#Jaccard figure
Jaccard_divergence <-
  data.frame(qute, jaccardvect1983, jaccardvect1985) %>%
  pivot_longer(2:3, names_to = "year", values_to = "differentation") %>%
  ggplot(aes(
    x = qute,
    y = differentation,
    group = year,
    color = year
  )) + geom_line() + geom_point()
Jaccard_divergence


