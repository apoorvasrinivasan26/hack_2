

library(tidyverse)
library(sqldf)
library(nlme)
library(lattice)
library(ggplot2)
library(haven)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend)
library(e1071)
df = read.csv("cleaneddata4.csv") %>%
  janitor::clean_names() %>%
  mutate(day = as.factor(day),
         dose = as.factor(dose),
         study = as.factor(new_seq_id))

#df = df %>%
 # filter(day == 0)
  

df_updated = read_csv("hack2_updated.csv")
#length(df$new_seq_id[which(df$study == 1 & df$day == 1 & df$dose ==0)])
#30480


#length(unique(df$new_seq_id[which(df$study == 1 & df$day == 1 & df$dose ==0)]))
##5080 proteins

#unique(df$dose[which(df$study == 1 & df$day == 0)])

#length(df$rfu2[which(df$study == 1 & df$day == 1 & df$dose ==0 & df$animal_new_id ==28)])

##spaghetti plot
df %>%
  filter(new_seq_id == 2) %>%
  ggplot(aes(x = day, y = percent_bwt_chag, group = animal_new_id)) +
  geom_line() +
  theme_bw()



  

compsym = gls(percent_bwt_chag ~ dose + day + day*dose + new_seq_id,
              data = df,
              correlation = corCompSymm(form = ~1 | day),
              method="REML", na.action = na.omit, control = list(singular.ok = TRUE))


#% bwt change is more at dose 0.1 day 3 vs 0.1 day 4



unique(df$study)


#df_wide = df %>%
#  spread(key = new_seq_id, value = rfu2)

write.csv(df_wide, "df_wide.csv")
 
#df_wide = dcast(df, animal_new_id + study ~ new_seq_id, value.var = "rfu2")


###taylors df-- linear regression

hack2_wide = read_csv("bslRFU_d7wht.csv")

hack2_wide = hack2_wide %>%
  filter(dose == 0.01) 

model = lm(day7_pchng_bwt ~ bsl_log10RFU, data = hack2_wide)
summary(model)


##PCA

wide_df = read_sas("final_wide.sas7bdat") 

pca_df = wide_df%>%
  select(-c(percent_bwt_chag)) 

pca_df[] <- lapply(pca_df, as.numeric)

  
pca= prcomp(pca_df, scale. = TRUE)

summary(pca)

pca$sdev

plot(pca)


##variance
plot(pca, type = "line")

summary(pca)
library(ggbiplot)

ggbiplot(pca, ellipse = TRUE)

# variance
pr_var = ( pca$sdev )^2 

# % of variance
prop_varex = pr_var / sum( pr_var )

# Plot
plot( prop_varex, xlab = "Principal Component", 
      ylab = "Proportion of Variance Explained", type = "b" )


# Scree Plot
plot( cumsum( prop_varex ), xlab = "Principal Component", 
      ylab = "Cumulative Proportion of Variance Explained", type = "b" )

##extracting pca

wide_df = cbind(wide_df, pca$x[, 1:20])

ggplot(pca_df, aes(PC1, PC2, col = protein_1620, fill = protein_1620)) +
  stat_ellipse(geom= "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")

adonis(pca_df[,1:] ~ iid)

###clustering



unsuper_df = wide_df %>%
  mutate(time = as.factor(time),
         dose =as.factor(dose),
         study = as.factor(study)) %>%
  select(-percent_bwt_chag)

bwt_chag = wide_df %>%
  select(percent_bwt_chag) %>%
  mutate(percent_bwt_chag = ifelse(percent_bwt_chag <=0, 1, 0))


#train_df = model_df[1:150,]
#test_df =model_df[151:228, ]


##hirarirchal model

data.scaled <- scale(unsuper_df)

# Calculate the (Euclidean) distances: data.dist
data.dist <- dist(data.scaled)

# Create a hierarchical clustering model: wisc.hclust
train.hclust <- hclust(data.dist, method = "complete")



# Cut tree so that it has 4 clusters: wisc.hclust.clusters
wisc.hclust.clusters <- cutree(train.hclust, k = 4)

# Compare cluster membership to actual diagnoses
hclust(wisc.hclust.clusters, bwt_chag)

plot(train.hclust)

summary(train.hclust)


##kmeans clustering

wisc.km <- kmeans(scale(unsuper_df), centers = 2, nstart = 20)



###PCA on random forest proteins

pca_rf = read_csv("md.top.baseline.csv") %>%
  select(-Variable) 

pca_rf$MD = as.numeric(pca_rf$MD)

pca_rf$VIMP = as.numeric(pca_rf$VIMP)

pca_new = prcomp(pca_df, scale. = TRUE) 

plot(pca_new, type = "line")

pr_var = ( pca_new$sdev )^2 

# % of variance
prop_varex = pr_var / sum( pr_var )

# Plot
plot( prop_varex, xlab = "Principal Component", 
      ylab = "Proportion of Variance Explained", type = "b" )


# Scree Plot
plot( cumsum( prop_varex ), xlab = "Principal Component", 
      ylab = "Cumulative Proportion of Variance Explained", type = "b")


pca_new = cbind(pca_rf, pca_new$x[, 1:20])




###k means clustering
set.seed(20)

wide_df = read_sas("final_wide.sas7bdat") %>%
  filter(time == 0)


final_df <- t(wide_df)

final_df = scale(final_df)

final_df = na.omit(final_df)
heatmap(final_df)

kmeans_df = final_df[6:5085, ]

kmean6 <- kmeans(kmeans_df, centers = 6, nstart = 25)
str(k2)

##use pearson distance


print(k2)

##vizuialization
fviz_cluster(kmean6, data = kmeans_df)

k_cluster = kmean6$cluster

k_cluster = data.frame(k_cluster)

write.csv(k_cluster, "kcluster_dat")

# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(kmeans_df, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")


##optimal k

set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(kmeans_df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


##there are 6 optimal clusters

##taylor's df
t_df = read_sas("animals_wide.sas7bdat") %>%
  select(-new_seqID)

kmeans_model = kmeans(t_df, centers = 6, nstart = 25)

fviz_cluster(kmeans_model, data = t_df)

prcomp(t_df, scale = TRUE)



##merging the two datasets

bs_day7 = read_csv("bslRFU_d7wht (1).csv") 


k_cluster = read_csv("kcluster_dat") 

bs_day7_kclust = sqldf("SELECT *
                       FROM bs_day7 as bs
                       LEFT JOIN k_cluster as k
                       ON bs.new_seqID = k.new_seq_id") %>%
  mutate(cluster = as.factor(cluster),
         day7_pchng_bwt = ifelse(day7_pchng_bwt <= -3, 1, ifelse((day7_pchng_bwt <= 0), 2, 3)),
         day7_pchng_bwt = as.factor(day7_pchng_bwt)) %>%
  select(-new_seqID) 

bs_day7_kclust$cluster = as.factor(bs_day7_kclust$cluster)

to_merge_df = bs_day7_kclust %>%
  select(cluster, iid, new_seq_id)
  
  
##pday7_pchng_bwt ==1 if less than -3, 2 if less than 0 3 if  0 and 4 if greater

clust_model = lm(day7_pchng_bwt ~ cluster + dose + study, data = bs_day7_kclust)
summary(clust_model)


range(bs_day7_kclust$day7_pchng_bwt[which(bs_day7_kclust$cluster == 1)])
# -6.5  2.9


prop_df = bs_day7_kclust %>%
  select(day7_pchng_bwt, cluster) 

prop.table(table(prop_df), 1)


###cluster 4 explains most variability

##picking proteins
c4_df = bs_day7_kclust %>%
    filter(cluster == 4) %>%
   filter(day7_pchng_bwt == 1) 


  
###analysis on clusters



##merging data for random forest

set.seed(123)
rf_merged_df = merge(wide_df[, 1:150], to_merge_df, on = new_seq_id)  %>%
  select(cluster, starts_with("protein_"))

library(randomForestSRC)
library(ramdomForest)
  
rf.model <-  rfsrc(cluster ~., data=rf_merged_df, ntree = 100, tree.err = TRUE, 
                   importance = TRUE, statistics = TRUE, predicted = TRUE, predicted.obb = TRUE, na.action='na.omit')


library (randomForest)
rf.model <- randomForest(cluster ~., data=rf_merged_df, ntree = 50, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)


# variable importance

# minimal depth
rf.model.minimal.depth <- var.select(cluster~.,data=rf_merged_df,object=rf.model)
md <- gg_minimal_depth(rf.model)

md.top.2 <- data.frame(Variable = md$varselect$names, MD = md$varselect$depth, VIMP =  md$varselect$vimp)
plot.gg_minimal_vimp(md)

###correlation plot
 
library(corrplot)

corr_df = wide_df %>%
  filter(dose == 0.4) %>%
  select(starts_with("protein_"))

corrplot(cor(corr_df[, 1000:1100]), method = "color", tl.cex = 0.5)



##do not use dependent variables in pca

##building model with pca_df

pca_df = wide_df %>%
  select(-starts_with("protein"))

colnames(pca_df)

##calculating the diff in protein exp at time 4 and 1
wide_df$protein_1[which(wide_df$time == 4)] - wide_df$protein_1[which(wide_df$time == 0)]


##hclust

hclust_df = wide_df %>%
  select(-percent_bwt_chag) 

hclust_df = as.data.frame(scale(hclust_df))


dist_mat <- dist(hclust_df, method = 'euclidean')

hclust_centroid <- hclust(dist_mat, method = 'centroid')
plot(hclust_centroid)


agnes = agnes(dist_mat, method = 'centroid')

cluster_cut = cutree(hclust_centroid, 3)

table(cluster_cut, wide_df$percent_bwt_chag)


##agnes
library(cluster)

t_df = as.data.frame(scale(t_df))

t_df = t_df[1:150, ]

dist_mat <- dist(t_df, method = 'euclidean')

hclust_ward <- hclust(dist_mat, method = 'average')

plot(hclust_ward)


  
###spliting data into training and test
  

set.seed(123) 
# Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(wide_df), size = floor(.75*nrow(wide_df)), replace = F)
train <- wide_df[sample, ]
test  <- wide_df[-sample, ]

##one way anova

bs_day7_kclust$cluster = as.factor(bs_day7_kclust$cluster)


df_cluster = bs_day7 %>%
  filter(study == 1) %>%
group_by(iid) %>%
  summarise(
    count = n(),
    mean = mean(day7_pchng_bwt, na.rm = TRUE),
    sd = sd(day7_pchng_bwt, na.rm = TRUE)
  )  %>%
  arrange(mean)

###5 out of 21 animals put on weight 


library("ggpubr")
ggboxplot(bs_day7_kclust, x = "iid", y = "day7_pchng_bwt", 
          color = "new_seq_id",
        #  order = c("1", "2", "3", "4", "5", "6"),
          ylab = "Weight", xlab = "iid")


##mean bwt change by new id-- not informative-- they all have the same mean 

protein_bwt_df  = bs_day7 %>%
  filter(study == 1) %>%
  group_by(new_seqID) %>%
    summarise(
      count = n(),
      mean = mean(day7_pchng_bwt, na.rm = TRUE),
      sd = sd(day7_pchng_bwt, na.rm = TRUE)
    )  %>%
    arrange(mean)


####vv informative on baseline rfu and seq id-- arranges proteins based on decs rfu values at baseline
bs_rfu_seq_id = bs_day7 %>%
  filter(study == 1) %>%
  group_by(new_seqID) %>%
    summarise(
      count = n(),
      mean = mean(bsl_log10RFU, na.rm = TRUE),
      sd = sd(bsl_log10RFU, na.rm = TRUE)
    )  %>%
    arrange(desc(mean))

##exporting df
write.csv(bs_rfu_seq_id, "bs_rfu_seq_id.csv")


bs_day7 = bs_day7 %>%
  mutate(pchng_bwt = ifelse(day7_pchng_bwt <= -3, 1, ifelse((day7_pchng_bwt <= 0), 2, 3)),
         pchng_bwt = as.factor(pchng_bwt)) 

res.aov <- aov(bsl_log10RFU ~ pchng_bwt, data = bs_day7)
# Summary of the analysis
summary(res.aov)


##found no significant difference in exp at baseline rfu and percent_bwt_change

TukeyHSD(res.aov)


"            diff          lwr          upr     p adj
2-1  0.002003154 -0.004191734 0.0081980411 0.7288653
3-1 -0.003124412 -0.009843711 0.0035948866 0.5203917
3-2 -0.005127566 -0.011090352 0.0008352204 0.1084779"


bs_day7$new_seqID = as.factor(bs_day7$new_seqID)

aov.seqid <- aov(day7_pchng_bwt ~ new_seqID, data = bs_day7)
# Summary of the analysis
summary(res.aov)

####following a paper- https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061002#ack


###machine learning--- final!!!!!!!

##SVM RFE

##from sigFeature package

##practice
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("sigFeature")

library(sigFeature)
#library(parallel)
library(SummarizedExperiment)

data(ExampleRawData, package="sigFeature") 
ExampleRawData



x <- t(assays(ExampleRawData)$counts)
y <- colData(ExampleRawData)$sampleLabels

pvals <- sigFeaturePvalue(x,y) 
hist(unlist(pvals),breaks=seq(0,0.08,0.0015),col="skyblue",
                                    xlab="p value",ylab="Frequency",main="")
##real

wide_df = read_sas("final_wide.sas7bdat") 
svm_df = wide_df %>%
  filter(dose == 0.01) %>%
  filter(!(time == 0)) %>%
  mutate(percent_bwt_chag = ifelse(percent_bwt_chag <= -1, 1, 0),
         percent_bwt_chag = as.factor(percent_bwt_chag)) %>%
  select(-study, -time, -dose, -iid, -log10_PK) %>%
  select(percent_bwt_chag, everything())


y_df = svm_df %>%
  select(percent_bwt_chag) %>%
  mutate(percent_bwt_chag = as.factor(percent_bwt_chag))
###removing subject 1026, 1017 cuz they put on weight


y_matrix = as.matrix(y_df) 
y_matrix =   as.numeric(y_matrix)

names(y_matrix) = c(1:18)

x_df = svm_df %>%
  select(-percent_bwt_chag)

x_matrix = as.matrix(x_df)

pvals <- sigFeaturePvalue(x_matrix, y_matrix)

hist(unlist(pvals),col="skyblue",breaks=seq(0, 1, 0.01),
                                    xlab="p value",ylab="Frequency", main="")


system.time(sigfeatureRankedList <- sigFeature(x_matrix, y_matrix))


write.table(sigfeatureRankedList, "sigfeatureRankedList_dose0.1_wologconc.txt")

n = nrow(x_matrix)

x_train = x_matrix[1:round(0.8*n),]
x_test = x_matrix[(31:n), ]

y_train = y_matrix[1:30]
y_test = y_matrix[31:n]

library(e1071)
sigFeature.model=svm(x_train[ ,sigfeatureRankedList[1:1000]], y_train,
                     type="C-classification", kernel="linear") 
summary(sigFeature.model)

pred <- predict(sigFeature.model, x_test[ ,sigfeatureRankedList[1:1000]]) 
table(pred,y_test)

"    y_test
pred 0 1
0 5 2
1 0 2"

system.time(featureRankedList <- svmrfeFeatureRanking(x_matrix, y_matrix))

write.table(featureRankedList, "featureRankedList_dose0.1_woconc.txt")



RFE.model=svm(x_train[ ,featureRankedList[1:1000]], y_train, type="C-classification", kernel="linear")
summary(RFE.model)

pred <- predict(RFE.model, x_test[ ,featureRankedList[1:1000]]) 
table(pred,y_test)

##at dose 0.4
"   y_test
pred 0 1
0 0 0
1 1 3"
##TP= 3, TN = 0, FP = 1, FN = 0



"    y_test
pred 0 1
0 3 1
1 0 2"

## accuracy = 3/4 = 83%


##rfe model at dose 0.1
"    y_test
pred 0 1
   0 1 1
   1 2 2"


"> pred <- predict(sigFeature.model, x_test[ ,sigfeatureRankedList[1:1000]]) 
> table(pred,y_test)
y_test
pred 0 1
0 3 1
1 0 2"


library(ROCR)
plot.roc(svm_df$percent_bwt_chag, pred)

pvalsigFe <- sigFeaturePvalue(x_matrix, y_matrix, 100, sigfeatureRankedList)
pvalRFE <- sigFeaturePvalue(x_matrix, y_matrix, 100, featureRankedList) 
par(mfrow=c(1,2))
hist(unlist(pvalsigFe),breaks=50, col="skyblue", main=paste("sigFeature"),
     xlab="p value") 

hist(unlist(pvalRFE),breaks=50, col="skyblue",
                          main=paste("SVM-RFE"), xlab="p value")

#In this section a comparison of p values produced by using top 100 features (using “sigFeature” and “SVM-RFE” algorithm) is shown by using scatter box plot.

mytitle<-'Box Plot'
boxplot(unlist(pvalsigFe), unlist(pvalRFE), main=mytitle,
        names=c("sigFeature", "SVM-RFE"),
        ylab="p value", ylim=c(min(unlist(pvalsigFe)), max(unlist(pvalRFE))))
stripchart(unlist(pvalsigFe), vertical=TRUE, method="jitter", add=TRUE, pch=16, col=c('green'))
stripchart(unlist(pvalRFE), vertical=TRUE, at=2, method="jitter", add=TRUE,
           pch=16, col=c('blue'))
grid(nx=NULL, ny=NULL, col="black", lty="dotted")

library("pheatmap")
library("RColorBrewer")


pheatmap(x_matrix[ ,featureRankedList[1:20]], scale="row", clustering_distance_rows="correlation")

pheatmap(x_matrix[ ,sigfeatureRankedList[1:20]], scale="row",clustering_distance_rows="correlation")

set.seed(1234)
results = sigFeature.enfold(x_matrix[, featureRankedList[1:1000]], y_matrix, "kfold" , 10)
data("results")
str(results[1])


FeatureBasedonFrequency <- sigFeatureFrequency(x_matrix[, featureRankedList[1:1000]], results, 50, 50, pf=FALSE) 
str(FeatureBasedonFrequency[1])

inputdata <- data.frame(y= y_df ,x=x_matrix[, featureRankedList[1:1000]])


sigCV = sigCVError(1:10, results , inputdata)

WritesigFeature(results, x_matrix[, featureRankedList[1:1000]])


##study pred at differnt doses from sigFeaturelist

#for dose 0.1

"> pred <- predict(sigFeature.model, x_test[ ,sigfeatureRankedList[1:1000]]) 
> table(pred,y_test)
y_test
pred 0 1
0 3 1
1 0 2"


##at dose 0.4
"   y_test
pred 0 1
0 0 0
1 1 3"
##TP= 3, TN = 0, FP = 1, FN = 0

##at dose 0.033



###same thing for different dose

wide_df = read_sas("final_wide.sas7bdat") 
svm_df = wide_df %>%
  filter(dose == 0.1) %>%
  filter(!(time == 0)) %>%
  mutate(percent_bwt_chag = ifelse(percent_bwt_chag <= -1, 1, 0),
         percent_bwt_chag = as.factor(percent_bwt_chag)) %>%
  select(-study, -time, -dose, -iid, -log10_PK) %>%
  select(percent_bwt_chag, everything())


y_df = svm_df %>%
  select(percent_bwt_chag) %>%
  mutate(percent_bwt_chag = as.factor(percent_bwt_chag))
###removing subject 1026, 1017 cuz they put on weight


y_matrix = as.matrix(y_df) 
y_matrix =   as.numeric(y_matrix)

names(y_matrix) = c(1:36)

x_df = svm_df %>%
  select(-percent_bwt_chag)

x_matrix = as.matrix(x_df)

pvals <- sigFeaturePvalue(x_matrix, y_matrix)

hist(unlist(pvals),col="skyblue",breaks=seq(0, 1, 0.01),
     xlab="p value",ylab="Frequency", main="")


system.time(sigfeatureRankedList <- sigFeature(x_matrix, y_matrix))


write.table(sigfeatureRankedList, "sigfeatureRankedList_dose0.0.033_wologconc.txt")

n = nrow(x_matrix)

x_train = x_matrix[1:round(0.8*n),]
x_test = x_matrix[30:n, ]

y_train = y_matrix[1:29]
y_test = y_matrix[30:n]

library(e1071)
sigFeature.model=svm(x_train[ ,sigfeatureRankedList[1:1000]], y_train,
                     type="C-classification", kernel="linear") 
summary(sigFeature.model)

pred <- predict(sigFeature.model, x_test[ ,sigfeatureRankedList[1:1000]]) 
table(pred,y_test)

"    y_test
pred 0 1
0 5 2
1 0 2"

system.time(featureRankedList <- svmrfeFeatureRanking(x_matrix, y_matrix))

write.table(featureRankedList, "featureRankedList_dose0.1_woconc.txt")



RFE.model=svm(x_train[ ,featureRankedList[1:1000]], y_train, type="C-classification", kernel="linear")
summary(RFE.model)

pred <- predict(RFE.model, x_test[ ,featureRankedList[1:1000]]) 
table(pred,y_test)




library(ROCR)
plot.roc(svm_df$percent_bwt_chag, pred)

pvalsigFe <- sigFeaturePvalue(x_matrix, y_matrix, 100, sigfeatureRankedList)
pvalRFE <- sigFeaturePvalue(x_matrix, y_matrix, 100, featureRankedList) 


write.table(pvalsigFe, "pvalsigFe.txt")
write.table(pvalRFE, "pvalRFE.txt")

par(mfrow=c(1,2))
hist(unlist(pvalsigFe),breaks=50, col="skyblue", main=paste("SVM-RFE and t-statistic"),
     xlab="p value") 

hist(unlist(pvalRFE),breaks=50, col="skyblue",
     main=paste("SVM-RFE"), xlab="p value")

#In this section a comparison of p values produced by using top 100 features (using “sigFeature” and “SVM-RFE” algorithm) is shown by using scatter box plot.

mytitle<-'Box Plot'
boxplot(unlist(pvalsigFe), unlist(pvalRFE), main=mytitle,
        names=c("SVM_RFE&tstat","SVM_RFE"),
        ylab="p value", ylim=c(min(unlist(pvalsigFe)), max(unlist(pvalRFE))))
stripchart(unlist(pvalsigFe), vertical=TRUE, method="jitter", add=TRUE, pch=16, col=c('green'))
stripchart(unlist(pvalRFE), vertical=TRUE, at=2, method="jitter", add=TRUE,
           pch=16, col=c('blue'))
grid(nx=NULL, ny=NULL, col="black", lty="dotted")

library("pheatmap")
library("RColorBrewer")


pheatmap(x_matrix[ ,featureRankedList[1:20]], scale="row", clustering_distance_rows="correlation")

pheatmap(x_matrix[ ,sigfeatureRankedList[1:20]], scale="row",clustering_distance_rows="correlation")

set.seed(1234)
results = sigFeature.enfold(x_matrix[, featureRankedList[1:1000]], y_matrix, "kfold" , 10)
data("results")
str(results[1])


FeatureBasedonFrequency <- sigFeatureFrequency(x_matrix[, featureRankedList[1:1000]], results, 50, 50, pf=FALSE) 
str(FeatureBasedonFrequency[1])

inputdata <- data.frame(y= y_df ,x=x_matrix[, featureRankedList[1:1000]])


sigCV = sigCVError(1:10, results , inputdata)

WritesigFeature(results, x_matrix[, featureRankedList[1:1000]])


library(corrplot)

top_proteins = svm_df %>%
  select(protein_2802, protein_1612, protein_1501, protein_1620, protein_3093, protein_2735, protein_2939, protein_2465, protein_4691,
         protein_4463, protein_4820, protein_4979, protein_5078, protein_2928, protein_493, protein_852)
corrplot(cor(top_proteins), method = "color", tl.cex = 0.5)

