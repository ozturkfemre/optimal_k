##########################
#### Required Packages ###
##########################
library(cluster)
library(NbClust)
library(readr)
library(readxl)
library(clValid)
library(factoextra)
library(tidyverse)
library(magrittr)
library(fpc)


##########################
### Required Functions ###
##########################

###############################
### function for silhouette ###
###############################

sigraph <- function(data) {
  k <- c(2:10)
  nb <- NbClust(data, min.nc = 2, max.nc = 10, index = "silhouette", method = "kmeans")
  s <- as.vector(nb$All.index)
  plot(k, s,xlab =  "Cluster number k",
       ylab = "Silhouette Score",
       main = "Silhouette Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(s==max(s)) + 1, lwd=1, col="red", lty="dashed")
}

#############################
### function for ch index ###
#############################

calingraph <- function(data) {
  k <- c(2:10)
  nb <- NbClust(data, min.nc = 2, max.nc = 10, index = "ch", method = "kmeans")
  ch <- as.vector(nb$All.index)
  plot(k, ch,xlab =  "Cluster number k",
       ylab = "Caliński - Harabasz Score",
       main = "Caliński - Harabasz Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(ch==max(ch)) + 1, lwd=1, col="red", lty="dashed")
}

#############################
### function for db index ###
#############################

dbgraph <- function(data) {
  k <- c(2:10)
  nb <- NbClust(data, min.nc = 2, max.nc = 10, index = "db", method = "kmeans")
  db <- as.vector(nb$All.index)
  plot(k, db,xlab =  "Cluster number k",
       ylab = "Davies-Bouldin Score",
       main = "Davies-Bouldin Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(db==min(db)) + 1, lwd=1, col="red", lty="dashed")
}


###############################
### function for dunn index ###
###############################

dunngraph <- function(data) {
  k <- c(2:10)
  dunnin <- c()
  for (i in 2:10) {
    dunnin[i] <- dunn(distance = dist(data), clusters = kmeans(data, i)$cluster)
  }
  dunnin <- dunnin[2:10]
  plot(k, dunnin, xlab =  "Cluster number k",
       ylab = "Dunn Index",
       main = "Dunn Plot", cex.main=1,
       col = "dodgerblue1", cex = 0.9 ,
       lty=1 , type="o" , lwd=1, pch=4,
       bty = "l",
       las = 1, cex.axis = 0.8, tcl  = -0.2)
  abline(v=which(dunnin==max(dunnin)) + 1, lwd=1, col="red", lty="dashed")
}

#############
### study ###
#############

# In this study, the most widely used methods for determining the number of clusters, namely Average Silhouette, Caliński - Harabasz, Davies - Bouldin and Dunn Index were used. The performances of these methods were compared by Rand Index and Meila's Variation of Information criteria on nine real data sets where the number of clusters was known in advance.

#########################
### k-means algorithm ###
#########################

# K-means is a popular unsupervised learning algorithm used for clustering. It partitions a set of n data points into k clusters, where k is a user-defined parameter. 
# The k-means algorithm minimizes the sum of squared distances between the data points and their closest cluster center to determine total within sum of squares for each cluster. 

#The algorithm works as follows:
# Initialization: Choose k initial cluster centers (centroids) randomly or using some other method, such as k-means++.

# Assignment: Assign each data point to the closest cluster center based on the Euclidean distance between the data point and the centroid.

# Recalculation: Calculate the mean of all the data points assigned to each cluster and use this mean as the new cluster center (centroid).

# Repeat steps 2 and 3 until the cluster assignments stop changing, or a maximum number of iterations is reached.

# Output: The final clustering solution is the k clusters, with each cluster being represented by its centroid.


# See [1], [2] for detailed information.

#########################
### Silhouette Method ###
#########################

# The silhouette method is a common evaluation metric used in unsupervised learning to measure the quality of a clustering solution. 
# The silhouette value for a particular data point is a measure of how similar that data point is to the other data points in its own cluster, compared to other clusters. 
# The silhouette value ranges from -1 to 1, where a value of 1 indicates a strong similarity to the other data points in its own cluster, and a value of -1 indicates a strong similarity to data points in another cluster.

# In essence, the silhouette method measures how well-defined a particular cluster is, and how well-separated it is from other clusters. 
# This makes it a useful tool for evaluating the effectiveness of different clustering algorithms or for determining the optimal number of clusters in a given dataset.

# For each observation the silhouette method can be calculated as follows:
# calculate cluster tightness: the average distance to all observations in the same cluster, which is called a(i).
# calculate cluster separation: calculate the minimum distance to all the data points in a different cluster, which is called b(i).
# calculate silhouette coefficient: calculate its silhouette value "s(i)" as the difference between the b(i) and a(i), divided by the maximum of these two distances: (b(i) - a(i)) / max(a(i), b(i))

# Finally, calculate the average silhouette value for all observations: mean(si)

# see [3] for detailed information

###########################
### Calinski - Harabasz ###
###########################

# The Calinski-Harabasz index (also known as the Variance Ratio Criterion) is a commonly used evaluation metric for comparing different clustering solutions in unsupervised learning.
# It is a ratio of the between-cluster variance and the within-cluster variance, and it is used to determine the number of clusters that should be used in a clustering solution.

# For each cluster calinski - harabasz method can be calculated as follows:
# calculate within cluster sum of squares WCSS: calculate the sum of the squared distances between each data point and its corresponding cluster center. 
# calculate between cluster sum of squares BCSS: calculate the sum of the squared distances between each cluster center and the overall mean of all the data points. 

# (BCSS / WCSS) - (n-k) / (k-1) k: total cluster number, n:total observation number

# The Calinski-Harabasz index ranges from 0 to infinity, with a higher value indicating a better clustering.

# see [4] for detailed information

########################
### Davies - Bouldin ###
########################

# The Davies-Bouldin index (DBI) is a measure of the similarity between the clusters in a clustering solution. 
# The DBI is calculated as the average similarity between each cluster and its most similar cluster, where the similarity between two clusters is defined as the maximum distance between any data point in one cluster and its closest data point in the other cluster.

# DBI is calculated as follows:
# calculate intra-cluster dispersion: calculate the mean distance of all the data points in each cluster to the cluster center.
# calculate separation criteria: calculate the Euclidean distance between the cluster centers for each pair of clusters.
# find the most similar cluster: for each pair of clusters, calculate the similarity d(i, j) between the two clusters, as defined in the previous answer.
# Sum up the maximum similarity between each cluster and its most similar cluster.
# Divide the sum by the number of clusters to obtain the Davies-Bouldin index

# The Davies-Bouldin index ranges from 0 to infinity, with a lower value indicating a better clustering solution. 
# A DBI of 0 indicates that there is no similarity between any two clusters, while a high DBI value indicates that there is a high level of similarity between some of the clusters. 

# For detailed information please see [5]

##################
### Dunn Index ###
##################

# The Dunn index is a measure of the compactness and separability of the clusters in a clustering solution. 
# It is calculated as the ratio of the minimum separation to the maximum diameter. 

# Dunn Index is calculated as follows:
# calculate minimum separation: calculate the smallest distance between the observations from two different clusters.
# calculate maximum diameter: calculate the maximum distance between the observations in the same cluster.

# The Dunn index ranges from 0 to infinity, with a higher value indicating a better clustering solution. 
# A value of 1 indicates that the clusters are perfectly separated and perfectly compact, while a low value indicates that the clusters are either not separated or not compact.

# for detailed information please see [6]

#######################
### Data Prepration ###
#######################

# The data preparation process for clustering is an important step to ensure that the clustering algorithm produces meaningful and accurate results. Some common steps involved in the data preparation process for clustering include:

# Data cleaning: This involves removing missing values, duplicates, or outliers from the data. This step helps to ensure that the clustering algorithm is not biased towards certain data points or clusters.

# Data standardization: This involves scaling the data so that all the features have a similar range and impact on the clustering results. 

# Dimensionality reduction: If the data has a large number of features, it may be helpful to reduce the dimensionality of the data using techniques like principal component analysis (PCA) or linear discriminant analysis (LDA). This can help to reduce noise and simplify the data, making it easier to interpret the results.

# Data transformation: This involves transforming the data into a different representation that may be more suitable for clustering. For example, transforming the data into a polar coordinate system or a logarithmic scale can help to highlight patterns or relationships that may not be apparent in the original data.

# Feature selection: This involves selecting a subset of the features that are most relevant for clustering. This can be done manually or using automated methods like mutual information or chi-squared tests.

# It's important to keep in mind that the data preparation process will vary depending on the specific requirements of the clustering problem. The most appropriate data preparation steps will depend on the type of data, the clustering algorithm being used, and the goals of the clustering analysis.

# please see [7] , [8] for detailed information

##################
### Rand Index ###
##################

# The Rand index is a measure of similarity between two clusterings, and can be used to compare a clustering result with a set of ground truth labels. 
# Given a clustering result with k clusters and a set of ground truth labels with c classes, the Rand index measures the proportion of pairs of data points that are either both assigned to the same cluster in the clustering result and to the same class in the ground truth labels, or both assigned to different clusters in the clustering result and to different classes in the ground truth labels.

# please see [9] for detailed information

########################################
### Meila's Variation of Information ###
########################################

# Meila's variation of information is a measure of the difference between two clusterings. 
# It is used to evaluate the quality of a clustering result compared to a set of ground truth labels. 
# The variation of information is calculated as the mutual information between the clustering result and the ground truth labels, normalized by the entropy of either the clustering result or the ground truth labels.

# please see [10] for detailed information

################################################################################
################################################################################
################################################################################
########################### wine dataset (3 clusters) ##########################
################################################################################
################################################################################
################################################################################


wine <- read_csv("~/Desktop/R/thesis/realdatasets/wine.data", 
                 col_names = FALSE) # data import
df <- wine[-1] # removing class variable

#################
### data prep ###
#################

corr <- cor(df, method = "pearson") # correlation
corrplot(corr, method = "color") # no significant correlation found, no need for pca

df <- scale(df) # the unit of measurement of each variable is different. data needs to be scaled

plot(df, col=wine$X1) # to see how clusters are separated.

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df)
# silhouette offered 3 clusters

################
### calinski ###
################

calingraph(df)
# ch offered 3 clusters

######################
### davies-bouldin ###
######################

dbgraph(df)
# db offered 3 clusters

##################
### dunn index ###
##################

dunngraph(df)
# dunn offered 3 clusters

################################
### Validation Measurements ####
################################

valids <- list()

##################
### silhouette ###
##################

wine_valid <- data.frame()
kmsidata <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)

df <- as.data.frame(df)
wine$X1 <- as.numeric(wine$X1)

wine_valid[1,1] <- "silhouette"
wine_valid[1,2] <- cluster.stats(d = dist(df), wine$X1 , kmsidata$cluster)$corrected.rand
wine_valid[1,3] <- cluster.stats(d = dist(df), wine$X1, kmsidata$cluster)$vi
##########
### ch ###
##########

kmchdata <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)

wine_valid[2,1] <- "ch"
wine_valid[2,2] <- cluster.stats(d = dist(df), wine$X1 , kmchdata$cluster)$corrected.rand
wine_valid[2,3] <- cluster.stats(d = dist(df), wine$X1, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)

wine_valid[3,1] <- "db"
wine_valid[3,2] <- cluster.stats(d = dist(df), wine$X1 , kmdbdata$cluster)$corrected.rand
wine_valid[3,3] <- cluster.stats(d = dist(df), wine$X1, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df, "kmeans", k = 3, nstart = 25, graph = F)

wine_valid[4,1] <- "dunn"
wine_valid[4,2] <- cluster.stats(d = dist(df), wine$X1 , kmdunndata$cluster)$corrected.rand
wine_valid[4,3] <- cluster.stats(d = dist(df), wine$X1, kmdunndata$cluster)$vi

names(wine_valid) <- c("method", "rand", "mvi")
valids$wine <- wine_valid

################################################################################
################################################################################
################################################################################                                          
########################## column3c dataset (3 clusters) #######################
################################################################################
################################################################################
################################################################################

column3c <- read_table("~/Desktop/R/thesis/realdatasets/column_3C.dat", 
                       col_names = FALSE) # data import

#################
### data prep ###
#################

names(column3c)
column3c %<>% mutate( X7 = case_when(
  X7 == "DH" ~ 1,
  X7 == "SL" ~ 2,
  X7 == "NO" ~ 3
) )
df1 <- column3c[-7]

corr <- cor(df1, method = "pearson")
corrplot(corr, method = "color") # significant correlation found, pca should be applied

data.pca <- prcomp(df1, center = TRUE, scale. = TRUE)
data.pca # 2 principle components are ok
df1 <- predict(data.pca)[,1:2]
df1 <- as.data.frame(df1)

plot(df1, col=column3c$X7) # to see how clusters are separated.

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df1)
# 2 clusters offered, it is not successful

################
### calinsky ###
################

calingraph(df1)
# 2 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df1)
# 9 clusters offered, it is not successful

##################
### dunn index ###
##################

dunngraph(df1)
# 10 clusters offered, it is not successful

################################
### Validation Measurements ####
################################

##################
### silhouette ###
##################

column_valid <- data.frame()
kmsidata <- eclust(df1, "kmeans", k = 2, nstart = 25, graph = F)

df1 <- as.data.frame(df1)
column3c$X7 <- as.numeric(column3c$X7)

column_valid[1,1] <- "silhouette"
column_valid[1,2] <- cluster.stats(d = dist(df1), column3c$X7, kmsidata$cluster)$corrected.rand
column_valid[1,3] <- cluster.stats(d = dist(df1), column3c$X7, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df1, "kmeans", k = 2, nstart = 25, graph = F)

column_valid[2,1] <- "ch"
column_valid[2,2] <- cluster.stats(d = dist(df1), column3c$X7, kmchdata$cluster)$corrected.rand
column_valid[2,3] <- cluster.stats(d = dist(df1), column3c$X7, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df1, "kmeans", k = 9, nstart = 25, graph = F)

column_valid[3,1] <- "db"
column_valid[3,2] <- cluster.stats(d = dist(df1), column3c$X7, kmdbdata$cluster)$corrected.rand
column_valid[3,3] <- cluster.stats(d = dist(df1), column3c$X7, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df1, "kmeans", k = 10, nstart = 25, graph = F)

column_valid[4,1] <- "dunn"
column_valid[4,2] <- cluster.stats(d = dist(df1), column3c$X7, kmdunndata$cluster)$corrected.rand
column_valid[4,3] <- cluster.stats(d = dist(df1), column3c$X7, kmdunndata$cluster)$vi

names(column_valid) <- c("method", "rand", "mvi")
valids$column3c <- column_valid

################################################################################
################################################################################
################################################################################                                          
############################   e.coli (8 clusters)   ###########################
################################################################################
################################################################################
################################################################################

ecoli <- read_table("~/Desktop/R/thesis/realdatasets/ecoli.data", 
                    col_names = FALSE) # data import
#################
### data prep ###
#################

df2 <- ecoli [-1]
df2 <- df2 [-8]

corr <- cor(df2, method = "pearson")
corrplot(corr, method = "color") # no significant correlation is found, no need for pca 

df2 <- scale(df2)

plot(df2, col= ecoli$site)

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df2) 
# 4 clusters offered, it is not successful

################
### calinsky ###
################

calingraph(df2)
# 4 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df2)
# 5 clusters offered, it is not successful

##################
### dunn index ###
##################

dunngraph(df2)
# 6 clusters offered, it is not successful

################################
### Validation Measurements ####
################################

colnames(ecoli) <- c("squence_name","mcg","gvh","lip","chg","aac","alm1","alm2","site")

ecoli %<>% mutate( site = case_when(
  site == "cp" ~ 1,
  site == "im" ~ 2,
  site == "imS" ~ 3,
  site == "imL" ~ 4,
  site == "imU" ~ 5,
  site == "om" ~ 6,
  site == "omL" ~ 7,
  site == "pp" ~ 8
) )

##################
### silhouette ###
##################

ecoli_valid <- data.frame()
kmsidata <- eclust(df2, "kmeans", k = 4, nstart = 25, graph = F)

ecoli_valid[1,1] <- "silhouette"
ecoli_valid[1,2] <- cluster.stats(d = dist(df2), ecoli$site, kmsidata$cluster)$corrected.rand
ecoli_valid[1,3] <- cluster.stats(d = dist(df2), ecoli$site, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df2, "kmeans", k = 4, nstart = 25, graph = F)

ecoli_valid[2,1] <- "ch"
ecoli_valid[2,2] <- cluster.stats(d = dist(df2), ecoli$site, kmchdata$cluster)$corrected.rand
ecoli_valid[2,3] <- cluster.stats(d = dist(df2), ecoli$site, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df2, "kmeans", k = 5, nstart = 25, graph = F)

ecoli_valid[3,1] <- "db"
ecoli_valid[3,2] <- cluster.stats(d = dist(df2), ecoli$site, kmdbdata$cluster)$corrected.rand
ecoli_valid[3,3] <- cluster.stats(d = dist(df2), ecoli$site, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df2, "kmeans", k = 6, nstart = 25, graph = F)

ecoli_valid[4,1] <- "dunn"
ecoli_valid[4,2] <- cluster.stats(d = dist(df2), ecoli$site, kmdunndata$cluster)$corrected.rand
ecoli_valid[4,3] <- cluster.stats(d = dist(df2), ecoli$site, kmdunndata$cluster)$vi

names(ecoli_valid) <- c("method", "rand", "mvi")
valids$ecoli <- ecoli_valid


################################################################################
################################################################################
################################################################################                                          
############################   iris (3 clusters)   #############################
################################################################################
################################################################################
################################################################################

iris <- get(data("iris", package = "datasets"))

df3 <- iris[1:4] # data import

#################
### data prep ###
#################

iris %<>% mutate(Species = case_when(
  Species == "setosa" ~ 1,
  Species == "versicolor" ~ 2,
  Species == "virginica" ~ 3
)) 

corr <- cor(df3, method = "pearson")
corrplot(corr, method = "color") # significant correlation found, pca should be applied

data.pca <- prcomp(df3, center = TRUE, scale. = TRUE)
summary(data.pca) # 2 pc is ok

df3 <- predict(data.pca)[,1:2]
df3 <- as.data.frame(df3)
plot(df3, col = iris$Species)

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df3) 
# 2 clusters offered, it is not successful

################
### calinsky ###
################

calingraph(df3)
# 9 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df3)
# 2 clusters offered, it is not successful

##################
### dunn index ###
##################

dunngraph(df3)
# 2 clusters offered, it is not successful

################################
### Validation Measurements ####
################################

df3 <- as.data.frame(df3)

##################
### silhouette ###
##################

iris_valid <- data.frame()
kmsidata <- eclust(df3, "kmeans", k = 2, nstart = 25, graph = F)

iris_valid[1,1] <- "silhouette"
iris_valid[1,2] <- cluster.stats(d = dist(df3), iris$Species, kmsidata$cluster)$corrected.rand
iris_valid[1,3] <- cluster.stats(d = dist(df3), iris$Species, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df3, "kmeans", k = 9, nstart = 25, graph = F)

iris_valid[2,1] <- "ch"
iris_valid[2,2] <- cluster.stats(d = dist(df3), iris$Species, kmchdata$cluster)$corrected.rand
iris_valid[2,3] <- cluster.stats(d = dist(df3), iris$Species, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df3, "kmeans", k = 2, nstart = 25, graph = F)

iris_valid[3,1] <- "db"
iris_valid[3,2] <- cluster.stats(d = dist(df3), iris$Species, kmdbdata$cluster)$corrected.rand
iris_valid[3,3] <- cluster.stats(d = dist(df3), iris$Species, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df3, "kmeans", k = 2, nstart = 25, graph = F)

iris_valid[4,1] <- "dunn"
iris_valid[4,2] <- cluster.stats(d = dist(df3), iris$Species, kmdunndata$cluster)$corrected.rand
iris_valid[4,3] <- cluster.stats(d = dist(df3), iris$Species, kmdunndata$cluster)$vi

names(iris_valid) <- c("method", "rand", "mvi")
valids$iris <- iris_valid


################################################################################
################################################################################
################################################################################                                          
############################   haberman(2 clusters)   ##########################
################################################################################
################################################################################
################################################################################


haberman <- read_csv("~/Desktop/R/thesis/realdatasets/haberman.data", 
                     col_names = FALSE)
df4 <- haberman[1:3]

#################
### data prep ###
#################

corr <- cor(df4, method = "pearson")
corrplot(corr, method = "color") # no correlation, no need for pca

df4 <- scale(df4)
plot(df4, col = haberman$X4)

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

silgraph(df4) 
# 4 clusters offered, it is not successful

################
### calinsky ###
################

calingraph(df4)
# 5 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df4)
# 5 clusters offered, it is not successful

############
### dunn ###
############

dunngraph(df4)
# 6 clusters offered, it is not successful

#################################
### Validation Measurements ####
################################

df4 <- as.data.frame(df4)

##################
### silhouette ###
##################

haberman_valid <- data.frame()
kmsidata <- eclust(df4, "kmeans", k = 4, nstart = 25, graph = F)

haberman_valid[1,1] <- "silhouette"
haberman_valid[1,2] <- cluster.stats(d = dist(df4), haberman$X4, kmsidata$cluster)$corrected.rand
haberman_valid[1,3] <- cluster.stats(d = dist(df4), haberman$X4, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df4, "kmeans", k = 5, nstart = 25, graph = F)

haberman_valid[2,1] <- "ch"
haberman_valid[2,2] <- cluster.stats(d = dist(df4), haberman$X4, kmchdata$cluster)$corrected.rand
haberman_valid[2,3] <- cluster.stats(d = dist(df4), haberman$X4, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df4, "kmeans", k = 5, nstart = 25, graph = F)

haberman_valid[3,1] <- "db"
haberman_valid[3,2] <- cluster.stats(d = dist(df4), haberman$X4, kmdbdata$cluster)$corrected.rand
haberman_valid[3,3] <- cluster.stats(d = dist(df4), haberman$X4, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df4, "kmeans", k = 6, nstart = 25, graph = F)

haberman_valid[4,1] <- "dunn"
haberman_valid[4,2] <- cluster.stats(d = dist(df4), haberman$X4, kmdunndata$cluster)$corrected.rand
haberman_valid[4,3] <- cluster.stats(d = dist(df4), haberman$X4, kmdunndata$cluster)$vi

names(haberman_valid) <- c("method", "rand", "mvi")
valids$haberman <- haberman_valid


################################################################################
################################################################################
################################################################################                                          
########################   wdbc dataset(2 clusters)   ##########################
################################################################################
################################################################################
################################################################################


wdbc <- read_csv("~/Desktop/R/thesis/realdatasets/wdbc.data", 
                 col_names = FALSE)

#################
### data prep ###
#################

wdbc <- wdbc[1:12]
df5 <- wdbc[3:12]

corr <- cor(df5, method = "pearson")
corrplot(corr, method = "color") # high correlation, pca needs to be applied

data.pca <- prcomp(df5, center = TRUE, scale. = TRUE) # 2 pc is ok
df5 <- predict(data.pca)[,1:2]

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df5) 
# 2 clusters offered, it is successful

################
### calinsky ###
################

calingraph(df5)
# 2 clusters offered, it is successful

######################
### davies-bouldin ###
######################

dbgraph(df5)
# 7 clusters offered, it is not successful

############
### dunn ###
############

dunngraph(df5)
# 6 clusters offered, it is not successful

#################################
### Validation Measurements ####
################################

df5 <- as.data.frame(df5)

names(wdbc) <- c("ID", "Diagnosis", "radius", "texture", "perimeter", "area", "smoothness", "compactness", "concavity", "concave.points", "symmetry", "fractal.dimension" )

wdbc %<>% mutate( Diagnosis = case_when(
  Diagnosis == "M" ~ 1,
  Diagnosis == "B" ~ 2
) )

##################
### silhouette ###
##################

wdbc_valid <- data.frame()
kmsidata <- eclust(df5, "kmeans", k = 2, nstart = 25, graph = F)

wdbc_valid[1,1] <- "silhouette"
wdbc_valid[1,2] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmsidata$cluster)$corrected.rand
wdbc_valid[1,3] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df5, "kmeans", k = 2, nstart = 25, graph = F)

wdbc_valid[2,1] <- "ch"
wdbc_valid[2,2] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmchdata$cluster)$corrected.rand
wdbc_valid[2,3] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df5, "kmeans", k = 7, nstart = 25, graph = F)

wdbc_valid[3,1] <- "db"
wdbc_valid[3,2] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmdbdata$cluster)$corrected.rand
wdbc_valid[3,3] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df5, "kmeans", k = 6, nstart = 25, graph = F)

wdbc_valid[4,1] <- "dunn"
wdbc_valid[4,2] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmdunndata$cluster)$corrected.rand
wdbc_valid[4,3] <- cluster.stats(d = dist(df5), wdbc$Diagnosis, kmdunndata$cluster)$vi

names(wdbc_valid) <- c("method", "rand", "mvi")
valids$wdbc <- wdbc_valid


################################################################################
################################################################################
################################################################################                                          
#################   breast tissue dataset (6 clusters)   #######################
################################################################################
################################################################################
################################################################################


breastissue <- read_excel("~/Desktop/R/thesis/realdatasets/breastissue.xlsx") # data import

#################
### data prep ###
#################

df6 <- breastissue[-1] 

breastissue %<>% mutate(Class = case_when(
  Class == "car" ~ 1,
  Class == "fad" ~ 2,
  Class == "mas" ~ 3,
  Class == "gla" ~ 4,
  Class == "con" ~ 5,
  Class == "adi" ~ 6
)) 

corr <- cor(df6, method = "pearson")
corrplot(corr, method = "color") # pca is necessary

data.pca <- prcomp(df6, center = TRUE, scale. = TRUE)
(data.pca$sdev)^2 # 2 pc

df6 <- predict(data.pca)[,1:2]

plot(df6, col = breastissue$Class)


#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df6) 
# 2 clusters offered, it is not successful

################
### calinski ###
################

calingraph(df6)
# 4 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df6)
# 4 clusters offered, it is not successful

#############
### dunn ####
#############

dunngraph(df6)
# 2 clusters offered, it is not successful

#################################
### Validation Measurements ####
################################

df6 <- as.data.frame(df6)

##################
### silhouette ###
##################

breast_valid <- data.frame()
kmsidata <- eclust(df6, "kmeans", k = 2, nstart = 25, graph = F)

breast_valid[1,1] <- "silhouette"
breast_valid[1,2] <- cluster.stats(d = dist(df6), breastissue$Class, kmsidata$cluster)$corrected.rand
breast_valid[1,3] <- cluster.stats(d = dist(df6), breastissue$Class, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df6, "kmeans", k = 4, nstart = 25, graph = F)

breast_valid[2,1] <- "ch"
breast_valid[2,2] <- cluster.stats(d = dist(df6), breastissue$Class, kmchdata$cluster)$corrected.rand
breast_valid[2,3] <- cluster.stats(d = dist(df6), breastissue$Class, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df6, "kmeans", k = 4, nstart = 25, graph = F)

breast_valid[3,1] <- "db"
breast_valid[3,2] <- cluster.stats(d = dist(df6), breastissue$Class, kmdbdata$cluster)$corrected.rand
breast_valid[3,3] <- cluster.stats(d = dist(df6), breastissue$Class, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df6, "kmeans", k = 2, nstart = 25, graph = F)

breast_valid[4,1] <- "dunn"
breast_valid[4,2] <- cluster.stats(d = dist(df6), breastissue$Class, kmdunndata$cluster)$corrected.rand
breast_valid[4,3] <- cluster.stats(d = dist(df6), breastissue$Class, kmdunndata$cluster)$vi

names(breast_valid) <- c("method", "rand", "mvi")
valids$breast_tissue <- breast_valid


################################################################################
################################################################################
################################################################################                                          
##################   appendicitis dataset (2 clusters)   #######################
################################################################################
################################################################################
################################################################################


appendicitis <- read_csv("~/Desktop/R/thesis/realdatasets/appendicitis.csv") # data import

df7 <- appendicitis[-8]

#################
### data prep ###
#################

corr <- cor(df7, method = "pearson")
corrplot(corr, method = "color")

data.pca <- prcomp(df7, center = TRUE, scale. = TRUE)
summary(data.pca)
(data.pca$sdev)^2 # 3 pc is ok


df7 <- predict(data.pca)[,1:3]


plot(appendicitis, col = appendicitis$Class)


#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df7) 
# 3 clusters offered, it is not successful

################
### calinski ###
################

calingraph(df7)
# 4 clusters offered, it is not successful

######################
### davies-bouldin ###
######################

dbgraph(df7)
# 10 clusters offered, it is not successful

############
### dunn ###
############

dunngraph(df7)
# 10 cluster offered, it is not successful

#################################
### Validation Measurements ####
################################

df7 <- as.data.frame(df7)

appendicitis %<>% mutate(Class = case_when(
  Class == 0 ~ 1,
  Class == 1 ~ 2
)) 

##################
### silhouette ###
##################

appen_valid <- data.frame()
kmsidata <- eclust(df7, "kmeans", k = 3, nstart = 25, graph = F)

appen_valid[1,1] <- "silhouette"
appen_valid[1,2] <- cluster.stats(d = dist(df7), appendicitis$Class, kmsidata$cluster)$corrected.rand
appen_valid[1,3] <- cluster.stats(d = dist(df7), appendicitis$Class, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df7, "kmeans", k = 4, nstart = 25, graph = F)

appen_valid[2,1] <- "ch"
appen_valid[2,2] <- cluster.stats(d = dist(df7), appendicitis$Class, kmchdata$cluster)$corrected.rand
appen_valid[2,3] <- cluster.stats(d = dist(df7), appendicitis$Class, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df7, "kmeans", k = 10, nstart = 25, graph = F)

appen_valid[3,1] <- "db"
appen_valid[3,2] <- cluster.stats(d = dist(df7), appendicitis$Class, kmdbdata$cluster)$corrected.rand
appen_valid[3,3] <- cluster.stats(d = dist(df7), appendicitis$Class, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df7, "kmeans", k = 10, nstart = 25, graph = F)

appen_valid[4,1] <- "dunn"
appen_valid[4,2] <- cluster.stats(d = dist(df7), appendicitis$Class, kmdunndata$cluster)$corrected.rand
appen_valid[4,3] <- cluster.stats(d = dist(df7), appendicitis$Class, kmdunndata$cluster)$vi

names(appen_valid) <- c("method", "rand", "mvi")
valids$appendicitis <- appen_valid


################################################################################
################################################################################
################################################################################                                          
######################   userknow dataset (2 clusters)   #########################
################################################################################
################################################################################
################################################################################


userknow <- userknow <- read_excel("~/Desktop/R/thesis/realdatasets/userknow.xlsx")
df8 <-  userknow[-6]

#################
### data prep ###
#################
unique(userknow$UNS)
userknow %<>% mutate(UNS =  case_when(
  UNS == "very_low" ~ 1,
  UNS == "High" ~ 2,
  UNS == "Low" ~ 3,
  UNS == "Middle" ~ 4
)) 

corr <- cor(df8, method = "pearson")
corrplot::corrplot(corr, method = "color")

df8 <- scale(df8)

plot(userknow, col = userknow$UNS)

#######################################
### Determination of Cluster Number ###
#######################################

##################
### silhouette ###
##################

sigraph(df8) 
# 9 clusters offered, it is not successful


################
### calinsky ###
################

calingraph(df8)
# 7 clusters offered, it is not successful


######################
### davies-bouldin ###
######################

dbgraph(df8)
# 9 clusters offered, it is not successful

############
### dunn ###
############

dunngraph(df8)
# 8 clusters offered it is not successful

#################################
### Validation Measurements ####
################################

df7 <- as.data.frame(df7)

##################
### silhouette ###
##################

userknow_valid <- data.frame()
kmsidata <- eclust(df8, "kmeans", k = 7, nstart = 25, graph = F)

userknow_valid[1,1] <- "silhouette"
userknow_valid[1,2] <- cluster.stats(d = dist(df8), userknow$UNS, kmsidata$cluster)$corrected.rand
userknow_valid[1,3] <- cluster.stats(d = dist(df8), userknow$UNS, kmsidata$cluster)$vi

##########
### ch ###
##########

kmchdata <- eclust(df8, "kmeans", k = 7, nstart = 25, graph = F)

userknow_valid[2,1] <- "ch"
userknow_valid[2,2] <- cluster.stats(d = dist(df8), userknow$UNS, kmchdata$cluster)$corrected.rand
userknow_valid[2,3] <- cluster.stats(d = dist(df8), userknow$UNS, kmchdata$cluster)$vi

##########
### db ###
##########

kmdbdata <- eclust(df8, "kmeans", k = 7, nstart = 25, graph = F)

userknow_valid[3,1] <- "db"
userknow_valid[3,2] <- cluster.stats(d = dist(df8), userknow$UNS, kmdbdata$cluster)$corrected.rand
userknow_valid[3,3] <- cluster.stats(d = dist(df8), userknow$UNS, kmdbdata$cluster)$vi

############
### dunn ###
############

kmdunndata <- eclust(df8, "kmeans", k = 8, nstart = 25, graph = F)

userknow_valid[4,1] <- "dunn"
userknow_valid[4,2] <- cluster.stats(d = dist(df8), userknow$UNS, kmdunndata$cluster)$corrected.rand
userknow_valid[4,3] <- cluster.stats(d = dist(df8), userknow$UNS, kmdunndata$cluster)$vi

names(userknow_valid) <- c("method", "rand", "mvi")
valids$userknow <- userknow_valid


                              ######################
##############################                      ############################
############################## Irıs and ecoli plots ############################
##############################                      ############################
                              ######################

library(ggplot2)

plot(df3, col = iris$Species, cex.main=1, cex = 0.9 ,
     pch=4,
     bty = "l",
     cex.axis = 0.8, tcl  = -0.2, main = "d5", xlab = "PC1", ylab="PC2")

typeof(df3)
df3 <- as.data.frame(df3)
iris$Species <- as.factor(iris$Species)
ggplot(df3, aes(x=PC1 , y=PC2)) + 
  geom_point(aes(col = iris$Species))+
  xlab("PC1") +
  ylab("PC2") +
  theme_minimal()+
  scale_color_manual(legend_title, values = c("green", "blue","red"))

ecoli$site <- as.factor(ecoli$site)

df2 <- as.data.frame(df2)

legend_title <- "class"
colorz <- c("skyblue1", "thistle1", "maroon1", "green", "goldenrod1", "moccasin", "purple4", "wheat")

ggplot(df2) + 
  geom_point(aes(x=mcg , y=gvh, color = ecoli$site))+
  xlab("mcg") +
  ylab("gvh") +
  theme_minimal()+ 
  scale_color_manual(legend_title, values =   colorz)


################################################################################
################################################################################
################################################################################

#### Averages #####

### silhouette

silrand <- c(valids$wine[1,2], valids$column3c[1,2], valids$ecoli[1,2], valids$iris[1,2],
             valids$haberman[1,2], valids$wdbc[1,2], valids$breast_tissue[1,2], valids$appendicitis[1,2], valids$userknow[1,2])

mean(silrand)


silmvi <-  c(valids$wine[1,3], valids$column3c[1,3], valids$ecoli[1,3], valids$iris[1,3],
             valids$haberman[1,3], valids$wdbc[1,3], valids$breast_tissue[1,3], valids$appendicitis[1,3], valids$userknow[1,3])

mean(silmvi)


### ch

chrand <- c(valids$wine[2,2], valids$column3c[2,2], valids$ecoli[2,2], valids$iris[2,2],
                       valids$haberman[2,2], valids$wdbc[2,2], valids$breast_tissue[2,2], valids$appendicitis[2,2], valids$userknow[2,2])

mean(chrand)


chmvi <-  c(valids$wine[2,3], valids$column3c[2,3], valids$ecoli[2,3], valids$iris[2,3],
             valids$haberman[2,3], valids$wdbc[2,3], valids$breast_tissue[2,3], valids$appendicitis[2,3], valids$userknow[2,3])

mean(chmvi)

### db

dbrand <- c(valids$wine[3,2], valids$column3c[3,2], valids$ecoli[3,2], valids$iris[3,2],
            valids$haberman[3,2], valids$wdbc[3,2], valids$breast_tissue[3,2], valids$appendicitis[3,2], valids$userknow[3,2])

mean(dbrand)


dbmvi <-  c(valids$wine[3,3], valids$column3c[3,3], valids$ecoli[3,3], valids$iris[3,3],
            valids$haberman[3,3], valids$wdbc[3,3], valids$breast_tissue[3,3], valids$appendicitis[3,3], valids$userknow[3,3])

mean(dbmvi)


### dunn

dunnrand <- c(valids$wine[4,2], valids$column3c[4,2], valids$ecoli[4,2], valids$iris[4,2],
            valids$haberman[4,2], valids$wdbc[4,2], valids$breast_tissue[4,2], valids$appendicitis[4,2], valids$userknow[4,2])

mean(dunnrand)


dunnmvi <-  c(valids$wine[4,3], valids$column3c[4,3], valids$ecoli[4,3], valids$iris[4,3],
            valids$haberman[4,3], valids$wdbc[4,3], valids$breast_tissue[4,3], valids$appendicitis[4,3], valids$userknow[4,3])

mean(dunnmvi)
