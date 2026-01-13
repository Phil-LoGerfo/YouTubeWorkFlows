


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

t_depmap_exp <- data.frame(t(depmap_exp))

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

# fix the colnames!
colnames(depmap_exp_set)  <- rownames(depmap_exp_yes)

gene_names_dep_map <- strsplit(rownames(depmap_exp_set), "..", fixed = TRUE)

gene_names_dep_map <- unlist(gene_names_dep_map)

genes_to_get <- seq(1, length(gene_names_dep_map), by = 2)

gene_names_dep_map <- gene_names_dep_map[genes_to_get]

rownames(depmap_exp_set) <- gene_names_dep_map

t_depmap_exp_set <- data.frame(t(depmap_exp_set))


dep_exp_model_data <- dep_model_data [intersect(colnames(depmap_exp_set), rownames(dep_model_data)),]  


#######################################
# Step 2: Read in expression and meta data from Expression Atlas database
# Go to website https://www.ebi.ac.uk/gxa/experiments/E-GEOD-100833/Downloads
# To download transcript expression data: Click on link "All Normalized Expressions for this experiment" 
# TO download meta data: Click on link "Experimental Design"
# Read using read.delim(directory/filename)
# Make SURE expression and meta data files sets are named correctly: marray_exp and marray_meta, respectively
###########################################################################


marray_exp <- read.delim("E-GEOD-100833-A-GEOD-13158-normalized-expressions.tsv")

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

marray_model_data <- read.delim("E-GEOD-100833-experiment-design.tsv")

dups <- grep("duplicate", marray_model_data$Assay)

marray_model_data_no_dups <- marray_model_data[-dups,]

marray_model_data <- marray_model_data_no_dups

rownames(marray_model_data) <- marray_model_data$Assay

marray_exp_model_data<- marray_model_data [intersect(colnames(marray_exp_set), rownames(marray_model_data)),]  

###########################################


dep_means <-  lapply(t_depmap_exp_set, mean)
dep_means_frame <- data.frame(unlist(dep_means))

dep_sdevs <-  lapply(t_depmap_exp_set, sd)
dep_sdevs_frame <- data.frame(unlist(dep_sdevs))

dep_mean_frame <- cbind(dep_means_frame, dep_sdevs_frame)
colnames(dep_mean_frame) <- c("mean", "sd")

dep_mean_frame <- dep_mean_frame[order(dep_mean_frame$mean),]
dep_mean_frame$cv <- dep_mean_frame$sd/dep_mean_frame$mean


dep_mean_frame$exp_standardized_to_sd <-  (dep_mean_frame$mean - min(dep_mean_frame$mean)) *
  (max(na.omit(dep_mean_frame$sd)) - 0) / (max(dep_mean_frame$mean) - min(dep_mean_frame$mean)) + 0

dep_mean_frame$exp_standardized_to_cv <-  (dep_mean_frame$mean - min(dep_mean_frame$mean)) *
  (max(na.omit(dep_mean_frame$cv)) - 0) / (max(dep_mean_frame$mean) - min(dep_mean_frame$mean)) + 0
View(dep_mean_frame)
dep_mean_frame <- dep_mean_frame[-c(1:5),]

write_tsv(dep_mean_frame, "data_files/dep_mean_frame.tsv")


dep_over_2 <- dep_mean_frame[which(dep_mean_frame$mean > 2),]
dep_over_1 <- dep_mean_frame[which(dep_mean_frame$mean > 1),]



#########################################################


marray_mean_exp <- lapply(t_marray_exp_set, mean)
marray_mean_exp <- data.frame(unlist(marray_mean_exp))
marray_mean_sd <- lapply(t_marray_exp_set, sd)
marray_mean_sd <- data.frame(unlist(marray_mean_sd))
marray_mean_frame <- cbind(marray_mean_exp, marray_mean_sd)
colnames(marray_mean_frame) <- c("mean", "sd")
marray_mean_frame <- marray_mean_frame[order(marray_mean_frame$mean),]
marray_mean_frame$cv <- marray_mean_frame$sd/marray_mean_frame$mean

marray_mean_frame$exp_standardized_to_sd <-  (marray_mean_frame$mean - min(marray_mean_frame$mean)) *
  (max(na.omit(marray_mean_frame$sd)) - 0) / (max(marray_mean_frame$mean) - min(marray_mean_frame$mean)) + 0

marray_mean_frame$exp_standardized_to_cv <-  (marray_mean_frame$mean - min(marray_mean_frame$mean)) *
  (max(na.omit(marray_mean_frame$cv)) - 0) / (max(marray_mean_frame$mean) - min(marray_mean_frame$mean)) + 0
View(marray_mean_frame)

write_tsv(marray_mean_frame, "data_files/marray_mean_frame.tsv")



######################################################################



plot(dep_mean_frame$mean)

plot(dep_mean_frame$sd)

plot(dep_mean_frame$cv)


plot(dep_mean_frame$mean, dep_mean_frame$sd)

loess_model <- loess(dep_mean_frame$sd ~ dep_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)



plot(marray_mean_frame$mean)

plot(marray_mean_frame$sd)

loess_model <- loess(marray_mean_frame$sd ~ c(1:length(marray_mean_frame$mean)))
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


plot(marray_mean_frame$cv)


plot(marray_mean_frame$mean, marray_mean_frame$sd)

loess_model <- loess(marray_mean_frame$sd ~ marray_mean_frame$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)





plot(marray_mean_frame$mean)

plot(marray_mean_frame$sd)


loess_model <- loess(marray_mean_frame$sd ~ c(1:length(marray_mean_frame$mean)))
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)

lines(marray_mean_frame$exp_standardized_to_sd, col = "red")
##############



plot(dep_mean_frame$mean)

plot(dep_mean_frame$sd)


loess_model <- loess(dep_mean_frame$sd ~ c(1:length(dep_mean_frame$mean)))
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)
lines(dep_mean_frame$exp_standardized_to_sd, col = "red")




plot(dep_mean_frame$sd)
plot(dep_mean_frame$cv)




points(gene_frame_qc_ordered_by_mean["GAPDH..2597.", "index"],
       gene_frame_qc_ordered_by_mean["GAPDH..2597.", "standardized_mean_to_point_3"], col = "blue", pch = 5)
















plot(dep_gene_summary$mean,dep_gene_summary$sd )
loess_model <- loess(dep_gene_summary$sd~ c(1:length(dep_gene_summary$mean)))
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)
loess_model$x[1]
loess_model$x[4]
loess_model <- loess(dep_gene_summary$sd~ length(dep_gene_summary$mean))
loess_model <- loess(dep_gene_summary$sd~ dep_gene_summary$mean)
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)
plot(dep_gene_summary$mean,dep_gene_summary$sd )
lines(loess_model$x, loess_model$fitted,col = "green", lwd = 2)


#####################################################




