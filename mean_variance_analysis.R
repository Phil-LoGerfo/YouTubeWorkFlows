


######################## Part 1 #########################################
####################################### 
#######################################
# Step 1: Read in expression and meta data from Expression Atlas database
# Go to website https://www.ebi.ac.uk/gxa/experiments/E-GEOD-100833/Downloads
# To download transcript expression data: Click on link "All Normalized Expressions for this experiment" 
# TO download meta data: Click on link "Experimental Design"
# Read using read.delim(directory/filename)
# Make SURE expression and meta data files sets are named correctly: marray_exp and marray_model_data, respectively
###########################################################################

marray_exp <- read.delim("c:/Users/loger/Documents/R_Projects/RNA_seq/E-GEOD-100833-A-GEOD-13158-normalized-expressions.tsv")

View(table(marray_exp$Gene.Name))

empties <- which(marray_exp$Gene.Name == "")

marray_exp <- marray_exp[-empties,]

View(table(marray_exp$Gene.Name))

dgcr5 <- which(marray_exp$Gene.Name == "DGCR5")

marray_exp <- marray_exp[-15502,]

rownames(marray_exp) <- marray_exp$Gene.Name

marray_exp_set <- marray_exp[,c(4:1748)]

t_marray_exp_set <- data.frame(t(marray_exp_set))

##############################################

marray_model_data <- read.delim("c:/Users/loger/Documents/R_Projects/RNA_seq/E-GEOD-100833-experiment-design.tsv")

#dups <- grep("duplicate", marray_model_data$Assay)

#marray_model_data_no_dups <- marray_model_data[-dups,]

#marray_model_data <- marray_model_data_no_dups

rownames(marray_model_data) <- marray_model_data$Assay

#dups <- grep("duplicate", rownames(t_marray_exp_set))

#t_marray_exp_set_no_dups <- t_marray_exp_set[-dups,]

#marray_exp_set <- data.frame(t(t_marray_exp_set_no_dups))

#t_marray_exp_set <- t_marray_exp_set_no_dups

marray_exp_model_data<- marray_model_data [intersect(colnames(marray_exp_set), rownames(marray_model_data)),]  


############################################################





######################## Part 2 #########################################
####################################### 
# Step 1: Read in expression and meta data from DepMap database
# Go to website https://depmap.org/portal/data_page/?tab=allData
# Make sure to use file under Version: DepMap Public 25Q3 
# Under the list of Primary Files, locate and download 
# OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
# and
# Model.csv
# Read the files into your environment using read.csv(directory/filename, row.names = 1)
# Make SURE expression and meta data files sets are named correctly: depmap_exp and dep_model_data, respectively
# depmap_exp <-read.csv(directory/filename, row.names = 1) 
# dep_model_data <-read.csv(directory/filename, row.names = 1) 
# From there you should be able to follow along with the video.
###########################################################################

depmap_exp <- 
  read.csv("C:/Users/loger/Documents/DepMap/Y25Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv", row.names = 1)

dep_model_data <- read.csv("C:/Users/loger/Documents/DepMap/Y25Q3/Model.csv", row.names = 1)

dim(depmap_exp)

length(unique(depmap_exp$ModelID))

View(table(depmap_exp$ModelID))

View(depmap_exp[which(depmap_exp$ModelID == "ACH-000672"),])

View(table(depmap_exp$IsDefaultEntryForMC))
View(table(depmap_exp$IsDefaultEntryForModel))

depmap_exp_yes <- depmap_exp[which(depmap_exp$IsDefaultEntryForModel == "Yes"),]

dim(depmap_exp_yes)

rownames(depmap_exp_yes) <- depmap_exp_yes$ModelID

depmap_exp_set <- t(depmap_exp_yes[,c(6:length(colnames(depmap_exp_yes)))])

depmap_exp_set <- data.frame(depmap_exp_set)

# fix the colnames! The underscores "_" got replaced with dots "." when transposing the table. 
colnames(depmap_exp_set)  <- rownames(depmap_exp_yes)


#Keep the GeneCards Symbols and split off the NCBI gene number, demarcated by ".." 
gene_names_dep_map <- strsplit(rownames(depmap_exp_set), "..", fixed = TRUE)

gene_names_dep_map <- unlist(gene_names_dep_map)

genes_to_get <- seq(1, length(gene_names_dep_map), by = 2)

gene_names_dep_map <- gene_names_dep_map[genes_to_get]

rownames(depmap_exp_set) <- gene_names_dep_map

t_depmap_exp_set <- data.frame(t(depmap_exp_set))


dep_exp_model_data <- dep_model_data [intersect(colnames(depmap_exp_set), rownames(dep_model_data)),]  




####################### Part 3 ##############################
#
# Make data frames that contain mean expression and variation of expression
# For both expression sets.
# Columns include: 
# mean 
# standard deviation (sd) 
# variance (var)
# coefficient of variation (cv)
# mean standardized to sd range (mean_standardized_to_sd)
# mean standardized to cv range (mean_standardized_to_cv)
###############################################

################### marray

marray_mean_exp <- lapply(t_marray_exp_set, mean)
marray_mean_exp <- data.frame(unlist(marray_mean_exp))

marray_mean_sd <- lapply(t_marray_exp_set, sd)
marray_mean_sd <- data.frame(unlist(marray_mean_sd))

marray_mean_frame <- cbind(marray_mean_exp, marray_mean_sd)
colnames(marray_mean_frame) <- c("mean", "sd")


marray_var <-  lapply(t_marray_exp_set, var)
marray_var_frame <- data.frame(unlist(marray_var))


marray_mean_frame$var <- marray_var_frame$unlist.marray_var.

# Order the frame by gene expression levels from lowest to highest.
marray_mean_frame <- marray_mean_frame[order(marray_mean_frame$mean),]

# Make a new column for Coefficient of Variation
marray_mean_frame$cv <- marray_mean_frame$sd/marray_mean_frame$mean

# standardize the mean values to be within the range of the sd values 
# Zi = (Xi – min(X)) / (max(X) – min(X))
marray_mean_frame$mean_standardized_to_sd <-  (marray_mean_frame$mean - min(marray_mean_frame$mean)) *
  (max(na.omit(marray_mean_frame$sd)) ) / (max(marray_mean_frame$mean) - min(marray_mean_frame$mean)) 


# standardize the mean values to be within the range of the cv values 
# Zi = (Xi – min(X)) / (max(X) – min(X))
marray_mean_frame$mean_standardized_to_cv <-  (marray_mean_frame$mean - min(marray_mean_frame$mean)) *
  (max(na.omit(marray_mean_frame$cv)) ) / (max(marray_mean_frame$mean) - min(marray_mean_frame$mean)) 

View(marray_mean_frame)



################### depmap


dep_means <-  lapply(t_depmap_exp_set, mean)
dep_means_frame <- data.frame(unlist(dep_means))

dep_sdevs <-  lapply(t_depmap_exp_set, sd)
dep_sdevs_frame <- data.frame(unlist(dep_sdevs))

dep_mean_frame <- cbind(dep_means_frame, dep_sdevs_frame)
colnames(dep_mean_frame) <- c("mean", "sd")

dep_var <-  lapply(t_depmap_exp_set, var)
dep_var_frame <- data.frame(unlist(dep_var))

dep_mean_frame$var <- dep_var_frame$unlist.dep_var.


# Order the frame by gene expression levels from lowest to highest.
dep_mean_frame <- dep_mean_frame[order(dep_mean_frame$mean),]

# Make a new column for Coefficient of Variation
dep_mean_frame$cv <- dep_mean_frame$sd/dep_mean_frame$mean


# standardize the mean values to be within the range of the sd values 
# Zi = (Xi – min(X)) / (max(X) – min(X))
dep_mean_frame$mean_standardized_to_sd <-  (dep_mean_frame$mean - min(dep_mean_frame$mean)) *
  (max(na.omit(dep_mean_frame$sd))) / (max(dep_mean_frame$mean) - min(dep_mean_frame$mean)) 


# standardize the mean values to be within the range of the cv values 
# Zi = (Xi – min(X)) / (max(X) – min(X))
dep_mean_frame$mean_standardized_to_cv <-  (dep_mean_frame$mean - min(dep_mean_frame$mean)) *
  (max(na.omit(dep_mean_frame$cv))) / (max(dep_mean_frame$mean) - min(dep_mean_frame$mean)) 

#View(dep_mean_frame)
#dep_mean_frame <- dep_mean_frame[-c(1:5),]

######################################################






######################################################################
############################# Part 4 #####################################
# Visualize and interpret the mean and variation of gene expression.
# Compare the differences between the datasets. 
# Take a look at some genes that stand out in the data. 

###########################################
####################################################### marray

View(marray_mean_frame)

plot(marray_mean_frame$mean)
plot(dep_mean_frame$mean)


plot(marray_mean_frame$mean,marray_mean_frame$var)
loess_model <- loess(marray_mean_frame$var ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


plot(marray_mean_frame$mean,marray_mean_frame$sd)
loess_model <- loess(marray_mean_frame$sd ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


test_gene <- "ACTB"
points(marray_mean_frame[test_gene, "mean"],
       marray_mean_frame[test_gene, "sd"], pch = 9, col = "blue", cex = 3)


plot(marray_mean_frame$mean, marray_mean_frame$cv)
loess_model <- loess(marray_mean_frame$cv ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)

#########################################


plot(marray_mean_frame$mean)

plot(marray_mean_frame$sd)

loess_model <- loess(marray_mean_frame$sd ~ c(1:length(marray_mean_frame$mean)))
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)

lines(marray_mean_frame$mean_standardized_to_sd, col = "red", lwd = 2)

#################################################################



##########################################################
########################################################
#####################################################

plot(dep_mean_frame$mean)


plot(dep_mean_frame$mean, dep_mean_frame$var)
loess_model <- loess(dep_mean_frame$var ~ dep_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


plot(dep_mean_frame$mean, dep_mean_frame$sd)
loess_model <- loess(dep_mean_frame$sd ~ dep_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)

test_gene <- "ACTB"
points(dep_mean_frame[test_gene, "mean"],
       dep_mean_frame[test_gene, "sd"], pch = 9, col = "blue", cex = 3)

test_gene <- "GAPDH"
points(dep_mean_frame[test_gene, "mean"],
       dep_mean_frame[test_gene, "sd"], pch = 9, col = "red", cex = 3)

test_gene <- "SPARC"

points(dep_mean_frame[test_gene, "mean"],
       dep_mean_frame[test_gene, "sd"], pch = 9, col = "blue", cex = 3)

samps <-  rownames((t_depmap_exp_set[which(t_depmap_exp_set[,test_gene]== 0),]))
samp_data <-(dep_exp_model_data[samps,])
samp_data_exp <- cbind(samp_data, t_depmap_exp_set[which(t_depmap_exp_set[,test_gene]  == 0 ), test_gene ])

#################

plot(dep_mean_frame$mean, dep_mean_frame$cv)
loess_model <- loess(dep_mean_frame$cv ~ dep_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)

View(dep_mean_frame)

test_gene <- "FSHB"

points(dep_mean_frame[test_gene, "mean"],
       dep_mean_frame[test_gene, "cv"], pch = 9, col = "blue", cex = 3)

samps <-  rownames((t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0),]))
samp_data <-(dep_exp_model_data[samps,])
samp_data_exp <- cbind(samp_data, t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0), test_gene ])

View(samp_data_exp)

test_gene <- "IL3"

points(dep_mean_frame[test_gene, "mean"],
       dep_mean_frame[test_gene, "cv"], pch = 9, col = "blue", cex = 3)

samps <-  rownames((t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0),]))
samp_data <-(dep_exp_model_data[samps,])
samp_data_exp <- cbind(samp_data, t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0), test_gene ])


############################################################## 


plot(marray_mean_frame$mean, marray_mean_frame$sd)
loess_model <- loess(marray_mean_frame$sd ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


test_gene <- "FABP1"
points(marray_mean_frame[test_gene, "mean"],
       marray_mean_frame[test_gene, "sd"], pch = 9, col = "blue", cex = 3)


plot(marray_mean_frame$mean, marray_mean_frame$cv)
loess_model <- loess(marray_mean_frame$cv ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


test_gene <- "MTTP"
points(marray_mean_frame[test_gene, "mean"],
       marray_mean_frame[test_gene, "cv"], pch = 9, col = "blue", cex = 3)


View(table(marray_exp_model_data$Sample.Characteristic.organism.part.))


#########################################################
#####################
part_names <- names(table(marray_exp_model_data$Sample.Characteristic.organism.part.))
part_nums <- table(marray_exp_model_data$Sample.Characteristic.organism.part.)


marray_mean_frame_biopsy_type <- NULL

for (i in 1:length(part_names))
  
  
{
  temp_names <- rownames(marray_exp_model_data[which(marray_exp_model_data$Sample.Characteristic.organism.part. == part_names[i]),])
  
  t_temp_exp_set <- t_marray_exp_set[temp_names,]
  
  
  temp_mean_exp <- lapply(t_temp_exp_set, mean)
  temp_mean_exp <- data.frame(unlist(temp_mean_exp))
  
  temp_mean_sd <- lapply(t_temp_exp_set, sd)
  temp_mean_sd <- data.frame(unlist(temp_mean_sd))
  
  temp_mean_frame <- cbind(temp_mean_exp, temp_mean_sd)
  colnames(temp_mean_frame) <- c("mean", "sd")
  
  temp_var <-  lapply(t_temp_exp_set, var)
  temp_var_frame <- data.frame(unlist(temp_var))
  temp_mean_frame$var <- temp_var_frame$unlist.temp_var.
  
  
  # Order the frame by gene expression levels from lowest to highest.
  temp_mean_frame <- temp_mean_frame[order(temp_mean_frame$mean),]
  
  # Make a new column for Coefficient of Variation
  temp_mean_frame$cv <- temp_mean_frame$sd/temp_mean_frame$mean
  
  # standardize the mean values to be within the range of the sd values 
  # Zi = (Xi – min(X)) / (max(X) – min(X))
  temp_mean_frame$exp_standardized_to_sd <-  (temp_mean_frame$mean - min(temp_mean_frame$mean)) *
    (max(na.omit(temp_mean_frame$sd)) ) / (max(temp_mean_frame$mean) - min(temp_mean_frame$mean)) 
  
  
  # standardize the mean values to be within the range of the cv values 
  # Zi = (Xi – min(X)) / (max(X) – min(X))
  temp_mean_frame$exp_standardized_to_cv <-  (temp_mean_frame$mean - min(temp_mean_frame$mean)) *
    (max(na.omit(temp_mean_frame$cv)) ) / (max(temp_mean_frame$mean) - min(temp_mean_frame$mean)) 
  
  temp_mean_frame$biopsy <- part_names[i]
  
  temp_mean_frame$nums <- part_nums[i]
  
  marray_mean_frame_biopsy_type <- rbind(marray_mean_frame_biopsy_type, temp_mean_frame)
  
}
  
  #################################################
  #### Talk about IL12B
  
  View(marray_mean_frame_biopsy_type)
  
  
  #########################################
  
  test_gene <- "REG1B"
  
  samps <-  rownames((t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0),]))
  samp_data <-(dep_exp_model_data[samps,])
  samp_data_exp <- cbind(samp_data, t_depmap_exp_set[which(t_depmap_exp_set[,test_gene] > 0), test_gene ])
  
  
  
  ######################### END
  
  