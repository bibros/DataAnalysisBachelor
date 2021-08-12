# http://rstudio-pubs-static.s3.amazonaws.com/420777_0d5281b35cdf41c0b0bcb5d30b608a30.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/CytoML/inst/doc/flowjo_to_gatingset.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/flowWorkspace/inst/doc/flowWorkspace-Introduction.html

### Concept of the script
# It reads in FlowJo worksapces with flow cytometry data and displays the data as graphs
### Contents
# -1 Needed packages
# 1a. Reading in single panels with data processing, visualization and export
# 1b. Loadings all exported csv files from the single panels, summarizing the data and visualizing them
# 2a. Reading in all workspaces from SYTOX Blue stainings, summarizing the data and visualizing it


############### -1 Needed packages ##########################################################################

install.packages("xlsx")
install.packages("naniar")
install.packages("data.table")
install.packages("ggpubr")
install.packages("ggsci")
install.packages("svDialogs")

source("https://bioconductor.org/biocLite.R")
BiocManager::install("openCyto")
BiocManager::install("ggcyto")
BiocManager::install("CytoML")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.13')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

require(ggpubr)
require(ggsci)
require(xlsx)
library(naniar)
library(data.table)

library(openCyto)
library(ggcyto)
library(CytoML)
library(flowWorkspace)
library(svDialogs)

pkgs <- c("remedy", "dplyr", "rstudioapi", "cowplot", "gplots", "ggplot2", "styler",
          "stringr", "RColorBrewer", "readxl", "devtools", "tidyverse", "scales",
          "useful", "renv", "patchwork")

for(i in 1:length(pkgs)){
  if(!require(pkgs[i], character.only = T)){
    install.packages(pkgs[i])
    require(pkgs[i], character.only = T)
  }else{
    require(pkgs[i], character.only = T)
  }
}


############### 1a Get Flow data ############################################################################

# select the input folder which holds the workspace and the samples
sampleFolder <- selectDirectory(
  caption = "Select folder for working directory", label = "Select",
  path = getwd()
)
setwd(sampleFolder)

# Loading in workspace and samples
flowWorkspacePath <- paste0(sampleFolder, "/", list.files(path = sampleFolder, pattern = "*.wsp"))
  x <- list.dirs(paste0(sampleFolder, "/Unmixed"), full.names = T)[-1]
flowSamplePath <- x[!grepl("Reference Group" , x)]
GraphName <- sub('.*/', '', sampleFolder)
GraphName <- substr(GraphName, 1, 7) # (nchar(GraphName[1])-13))
stained <- 2

ws <- open_flowjo_xml(flowWorkspacePath)
print(ws)

gs <- flowjo_to_gatingset(ws,name = 1, path = flowSamplePath)  # takes all FCS files which are in the "all samples" group AND the Tumor folder
gs_pop_get_count_fast(gs) # <- these are the raw cell counts within the fcs files

#view(gs_pop_get_count_fast(gs))
#autoplot(gs[[stained]], "Cells/Single Cells/Live cells/Immune cells/")
#plot <- autoplot(gs[[stained]], c("Cells/Single Cells/Live cells/Immune cells", "Cells/Single Cells/Live cells/Fibro Endo", "Cells/Single Cells/Live cells/Tumor"))
#plot$`Cells`
#plot(gs)     # plot gating strategy




#-------------- Process Flow data for immune panel ----------------------------------------------------------

gs_pop_get_count_fast(gs[stained])$Count

# create new variable with essential data for calculating percentages
percentages <- data.frame(Population = rep(gs_pop_get_count_fast(gs[stained])$Population),
                          Count = rep(gs_pop_get_count_fast(gs[stained])$Count),
                          Parent = rep(gs_pop_get_count_fast(gs[stained])$ParentCount),
                          Percent = (rep(gs_pop_get_count_fast(gs[stained])$Count))/(rep(gs_pop_get_count_fast(gs[stained])$ParentCount))*100)
#head(percentages)

# get names of gated populations
pop <- gs_pop_get_count_fast(gs[stained])$Population
pop1 <- str_extract_all(pop, "[^/](.*?)+[^/]+", simplify = TRUE)
pop1[pop1 == ""] <- NA
#View(pop1)
lastValue <- function(pop1)   tail(pop1[!is.na(pop1)], 1)


# safe gated populations as new variable and adjust names
newGates <- apply(pop1, 1, lastValue)

newGates[19] = "Other 2"
newGates[26] = "Other 3"
newGates[33] = "Tc ctotoxic/exhausted"
newGates[34] = "Tc inflammatory"
newGates[37] = "T conv ctotoxic/exhausted"
newGates[38] = "T conv inflammatory"
#View(newGates)

percentages["Population"] <- newGates
percentages["PercTotal"] <- percentages$Count/percentages[3,2]*100   # % cells of live cells

plot <- percentages[4:39,]



#-------------- Plot and export data ------------------------------------------------------------------------

plot$Population <- factor(plot$Population, levels = plot$Population)

ggplot() +
  geom_col(data = plot, width = .2, aes(x = Population, y = PercTotal),
           color = "black") +
  
  ylab("Cell count (% of live cells)") +
  xlab("Cell types") +
  theme_classic() +
  
  ggtitle("TILs") +
  
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# esport as png and excel file
setwd("../_Export/Immune")
ggsave(plot = last_plot(), paste(GraphName, " Immune", ".png", sep = ""),
       width = 20, height = 10, units = "cm",
       bg = "transparent")
       # for date: paste(GraphName, format(Sys.time(), " %d-%m-%Y")

excelExport <- pop1
excelExport <- cbind(excelExport, percentages[,2:5])

write_excel_csv(excelExport,
                file = paste(GraphName, " Immune", ".csv", sep = ""))




############### 1b Summarize data ###########################################################################

# get all exported csv files from the immune panels
setwd("c:/Users/Bianca/Documents/Uni/6) Bachelorarbeit/6 Data/Flow Jo Analysis/Panel Stainings_Tumor_Immune/")
folder <- "_Export/Immune/"      # path to folder that holds the .csv files
file_list <- list.files(path = folder, pattern = "*.csv") # create list of all .csv files in folder
View(file_list)

#file_list <- file_list[c(11,12)]  # choose specific files

GateNames <- read_csv("_Export/GateNames_Immune2.csv")
View(GateNames)
samples <- GateNames

# each file is read, sample name is generated and all "percentages of total live cells" are combined in one data frame
for (i in 1:length(file_list)){
         sampleFile <- read.csv(paste(folder, file_list[i], sep=''))
         
         names(sampleFile)[1] = "X1"
         gatesOfSamples <- sampleFile %>% select(c(1:9))
         test <- function(gatesOfSamples)   tail(gatesOfSamples[!is.na(gatesOfSamples)], 1)
         sampleFile <- sampleFile %>% rename_at(13,~str_sub(file_list[i], 1, 7))
         
         samples <- cbind(samples, sampleFile[13])
         samples <- samples %>% rename_at(1, ~"Population")
}
#View(samples)


# selection of specific populations for plotting
SumPop1 <- samples[c(4,5,40),] # tumor/ fibro/ immune
        selection <- "Tumor-Fibro-Immune"
        height <- 8
SumPop1 <- samples[c(6,9,15,25,30,35),]
        selection <- "MajorClasses" # immune major cell classes
        height <- 15
SumPop1 <- samples[c(10,11,16:18,24,31,32,36,39),] # minor immune cell classes
        selection <- "MinorClasses"
        height <- 15
#-------------- 
SumPop1 <- samples[c(4,6,10:13,16:18,21,22,24,25,28,31,33,34,36,38:40),] # All end gates
        selection <- "AllEndPops"
        height <- 18
SumPop1 <- samples[c(4:40),]
        selection <- "AllGates"
        height <- 25

# transformation needed for plotting
SumPop <- as.data.table(SumPop1)
SumPop <- melt(SumPop, id.vars = c("Population"), variable.name = "TumorID", value.name = "PercOfTotal")

#
#### Plot all tumor samples as bars next to each other ----------------------

ggplot(data = SumPop, width = .2, aes(x = Population, y = PercOfTotal, fill = TumorID)) +
  geom_col(position = position_dodge2(preserve = "single"),
           color = "black") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm')) +
  labs(title = "Cells in tumors", y = "Percent of total live cells", x = "Cell types")

ggsave(plot = last_plot(), path = "_Export", paste("AllTumors_Immune_", selection, ".png", sep = ""),
       width = 20, height = height, units = "cm",
       bg = "transparent")



#-------------- Create means for entities with tumor samples as single values ----------

# create coloumn with corresponding entity for each tumor samples
SumPop[, Entity:=as.character(SumPop$TumorID)]
SumPop2 <- SumPop

for (i in 1:nrow(SumPop)){      # Get entities
  SumPop2[i,4] <- str_sub(SumPop[i,4], start = 0, end = 3)
}

SumPop3 <- aggregate(x = SumPop2$PercOfTotal,
           by = list(SumPop2$Entity, SumPop2$Population),
           FUN = mean) %>% rename_at(1:3, ~c("Entity", "Population", "MeanPerc"))
new <- SumPop2 %>%      # get sample numbers
  group_by(Entity) %>%
  summarise(no_samples = n_distinct(TumorID))
SumPop3 <- cbind(SumPop3, new[,2])



#-------------- Plot entity means, One graph -------------------------------------------

ggplot(data = SumPop3, width = .2, aes(x = Population, y = MeanPerc, fill = Entity)) +
  geom_col(position = position_dodge2(preserve = "single"),
           color = "black") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cells in tumors", y = "Percent of total live cells", x = "Cell types") +
  geom_point(data = SumPop2, aes(x = Population, y = PercOfTotal), #, shape = Entity), 
             color = "black", position=position_dodge(.9))  #+ # Punkte mit Einzelwerten
  #scale_shape_manual(values = c(8,15,16,17,18,19,9))


ggsave(plot = last_plot(), path = "_Export", paste("EntityMeans_OneGraph_Immune_", selection, ".png", sep = ""),
       width = 20, height = height, units = "cm",
       bg = "transparent")


#-------------- Plot entity means, Multiple graphs -------------------------------------

ggplot(data = SumPop3, width = .2, aes(x = Entity, y = MeanPerc, fill = Entity)) +
  facet_wrap(~Population, strip.position = "bottom", nrow = 2) +
  geom_col(position = position_dodge(width = 1),
           color = "black") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key.size = unit(0.5, 'cm'),
        text = element_text(size = 16),
        panel.background = element_rect(fill = "white")) +
  labs(title = " ", y = "Percent of total live cells", x = "Cell types") +
  geom_point(data = SumPop2, aes(x = Entity, y = PercOfTotal), #, shape = Entity), 
             color = "black", position=position_dodge(.9))  #+ # Punkte mit Einzelwerten
  #scale_shape_manual(values = c(8,15,16,17,18,19,9))

ggsave(plot = last_plot(), path = "_Export", paste("EntityMeans_MultipleGraphs_Immune_", selection, ".png", sep = ""),
       width = 22, height = height, units = "cm",
       bg = "transparent")

##




############### 2a Sytox stainings ##########################################################################

# get input directory
sampleFolder <- selectDirectory(
  caption = "Select folder for working directory", label = "Select",
  path = getwd()
)
setwd(sampleFolder)

# get all workspaces
file_list <- list.files(recursive = T, pattern = "*.wsp") # create list of all .csv files in folder
View(file_list)

files <- str_extract_all(file_list, "[^/](.*?)+[^/]+", simplify = TRUE)
files[files == ""] <- NA
colnames(files) = c("folder", "workspace")

# read in each workspaces, calculate percentages and export png/excel file for each sample
for (i in 1:length(file_list)){
  flowWorkspacePath <- file_list[i]
  
  x <- list.dirs(paste0(files[i,1], "/Unmixed"), full.names = T)[-1]
  flowSamplePath <- x[!grepl("Reference Group" , x)]
  
  GraphName <- str_sub(files[i, 2], 1, length(files[i, 2])-6)
  
  stained <- 1
  
  ws <- open_flowjo_xml(flowWorkspacePath)
  print(ws)
  
  gs <- flowjo_to_gatingset(ws,name = 1, path = flowSamplePath)  # takes all FCS files which are in the "all samples" group AND the Tumor folder
  gs_pop_get_count_fast(gs) # <- these are the raw cell counts within the fcs files
  
  percentages <- data.frame(Population = rep(gs_pop_get_count_fast(gs[stained])$Population),
                            Count = rep(gs_pop_get_count_fast(gs[stained])$Count),
                            Parent = rep(gs_pop_get_count_fast(gs[stained])$ParentCount),
                            Percent = (rep(gs_pop_get_count_fast(gs[stained])$Count))/(rep(gs_pop_get_count_fast(gs[stained])$ParentCount))*100)
  
  pop <- gs_pop_get_count_fast(gs[stained])$Population
  pop1 <- str_extract_all(pop, "[^/](.*?)+[^/]+", simplify = TRUE)
  pop1[pop1 == ""] <- NA
  
  lastValue <- function(pop1)   tail(pop1[!is.na(pop1)], 1)
  newGates <- apply(pop1, 1, lastValue)
  
  percentages["Population"] <- newGates
  percentages["PercTotal"] <- c(percentages$Count[1:4]/percentages[1,2]*100,  # % cells of total cells
                                percentages$Count[5:6]/percentages[4,2]*100)  # % of lymphocytes
  
  excelExport <- cbind(pop1, percentages[,2:5])
  write_excel_csv(excelExport, file = paste("_Export_csvFiles/", GraphName,
                                            ".csv", sep = ""))
}


folder <- "_Export_csvFiles/"      # path to folder that holds the .csv files
file_list <- list.files(path = folder, pattern = "*.csv") # create list of all .csv files in folder
#View(file_list)

GateNames <- read_csv("GateNames_Sytox.csv")
samples <- GateNames
samples <- samples %>% rename_at(1, ~"Population")

samples[3,] = "Live total cells"
samples[5,] = "Single cells lymphocytes"
samples[6,] = "Live lymphocytes"

# adjust names for each sample and combine percentages from each sample into one single data frame
for (i in 1:length(file_list)){
  sampleFile <- read.csv(paste(folder, file_list[i], sep=''))
  
  gatesOfSamples <- sampleFile %>% rename_at(1, ~"X1") %>% select(c(1:3))
  test <- function(gatesOfSamples)   tail(gatesOfSamples[!is.na(gatesOfSamples)], 1)
  sampleFile <- sampleFile %>% rename_at(7,~str_sub(file_list[i], 1, length(file_list[i])-6)) # changes colname for percentage of total cells
  
  samples <- cbind(samples, sampleFile[7])
}

# select populations of interest
SumPop1 <- samples[c(3),] # 3 total live cells / 4 lymphocytes / 6 live lyphocytes
SumPop <- as.data.table(SumPop1)
SumPop <- melt(SumPop, id.vars = c("Population"), variable.name = "TumorID", value.name = "PercOfTotal")


#-------------- Create means for entities with tumor samples as single values ----------

# create coloumn with entities
SumPop[, Entity:=as.character(SumPop$TumorID)]
SumPop2 <- SumPop

for (i in 1:nrow(SumPop)){      # Get entities
  SumPop2[i,4] <- str_sub(SumPop[i,4], start = 0, end = 3)
}

SumPop3 <- aggregate(x = SumPop2$PercOfTotal,
                     by = list(SumPop2$Entity, SumPop2$Population),
                     FUN = mean) %>% rename_at(1:3, ~c("Entity", "Population", "MeanPerc"))

# calculate standard deviation
sd_mean <- SumPop2 %>%
  group_by(Entity) %>%
  summarise( 
    n = n(),
    mean = mean(PercOfTotal),
    sd = sd(PercOfTotal)
  )
  

#-------------- Plot entity means, One graph -------------------------------------------

ggplot(data = SumPop3, width = .2, aes(x = Entity, y = MeanPerc, fill = Entity)) +
  geom_col(position = position_dodge2(preserve = "single"),
           color = "black") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        text = element_text(size = 14)) +
  labs(y = "Live cells \n [% of total cells]", x = "Entities") +
  geom_point(data = SumPop2, aes(x = Entity, y = PercOfTotal), size = 1.2, 
             color = "black", position=position_dodge(.9)) +  # Punkte mit Einzelwerten
  geom_errorbar(data = SumPop3, aes(x = Entity, ymin = sd_mean$mean-sd_mean$sd, ymax = sd_mean$mean+sd_mean$sd),
                width = 0.4, colour = "black", alpha = 0.9, size = 1)
  #title = "Live cells in tumors", 

title <- dlgInput("Enter the name of your graph for the export")$res

ggsave(plot = last_plot(), path = "_Export_graphs/", paste(title, ".png", sep = ""),
       width = 10, height = 7, units = "cm",
       bg = "transparent")




############### 3 Comparisons of immune panels / tumor vs. organoids ########################################

# get input directory
sampleFolder <- selectDirectory(
  caption = "Select folder for working directory", label = "Select",
  path = getwd()
)
setwd(sampleFolder)

folderTumor <- "Panel Stainings_Tumor/_Export/Immune/"      # path to folder that holds the .csv files
folderOrganoid <- "Panel Stainings_Organoids/_Export/"      # path to folder that holds the .csv files
file_listTumor <- list.files(path = folderTumor, pattern = "*.csv") # create list of all .csv files in folder
file_listOrgan <- list.files(path = folderOrganoid, pattern = "*.csv") # create list of all .csv files in folder

tumor <- NULL
organoid <- NULL

for (i in 1:length(file_listTumor)){
  tumor <- rbind(tumor, str_sub(file_listTumor[i], 1, 7))
}
for (i in 1:length(file_listOrgan)){
  organoid <- rbind(organoid, str_sub(file_listOrgan[i], 1, 7))
}

for (i in 1:length(file_listOrgan)){
  if (any(tumor == organoid[i]) == T){
    id <- organoid[i]
    fileTumor <- read.csv(paste(folderTumor, file_listTumor[which(tumor == organoid[i])], sep=''))
    fileOrgan <- read.csv(paste(folderOrganoid, file_listOrgan[i], sep=''))
    
    samples <- read_csv("Panel Stainings_Tumor/_Export/GateNames_Immune.csv")
    samples <- cbind(samples, fileTumor[10], fileOrgan[10])
    names(samples)[1:3] <- c("Population", "TumorCount", "OrganoidsCount")
    samples["Tumor"] <- samples$Tumor/samples[5,2]*100   # % cells of immune cells
    samples["Organoids"] <- samples$Organoids/samples[5,3]*100   # % cells of immune cells
    
    SumPop1 <- samples[c(6,9,15,20,24,25,30,35),c(1,4,5)] # select all immune cell gates
    SumPop <- as.data.table(SumPop1)
    SumPop <- melt(SumPop, id.vars = c("Population"), variable.name = "Tissue", value.name = "PercOfImmune")

#   View(SumPop)
    ggplot(data = SumPop, width = .2, aes(x = Tissue, y = PercOfImmune, fill = Population)) +
      geom_col(position = "stack", # position_dodge2(preserve = "single"),
               color = "black") +
      theme(panel.background = element_rect(fill = "white"),
            axis.line = element_line(color = "black"),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            text = element_text(size = 16),
            legend.key.size = unit(0.5, 'cm')) +
      labs(title = paste("Cells in", id, "tissues"), y = "Percent of total immune cells", x = "Cell types")
    
    setwd("c:/Users/Bianca/Documents/Uni/6) Bachelorarbeit/6 Data/Flow Jo Analysis/_Tumor-Organoid_comparison/")
    ggsave(plot = last_plot(), paste(id, "_stacked", ".png", sep = ""),
           width = 15, height = 25, units = "cm",
           bg = "transparent")

    ggplot(data = SumPop, width = .2, aes(x = Tissue, y = PercOfImmune, fill = Population)) +
      geom_col(position = "stack", # position_dodge2(preserve = "single"),
               color = "black") +
      theme(panel.background = element_rect(fill = "white"),
            axis.line = element_line(color = "black"),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text.y = element_text(angle = 90, hjust = 0.5),
            text = element_text(size = 18),
            legend.position = "bottom",
            legend.key.size = unit(0.5, 'cm')) +
      labs(title = paste("Cells in", id, "tissues"), y = "Percent of total immune cells", x = "Condition") +
      coord_flip()
    
    ggsave(plot = last_plot(), paste(id, "_stacked_quer", ".png", sep = ""),
           width = 25, height = 10, units = "cm",
           bg = "transparent")
  }
}





############### 4 Comparisons of immune panels / experiments ################################################

# get input directory
sampleFolder <- selectDirectory(
  caption = "Select folder for working directory", label = "Select",
  path = getwd()
)
setwd(sampleFolder)

folder <- selectDirectory(
  caption = "select folder with sample data", label = "select",
  path = getwd()
)

fileList <- list.files(path = folder, recursive = T, pattern = "*.csv")
View(fileList)

sample <- NULL

for (i in 1:length(fileList)){
  x <- sub('.*/', '', fileList[i])
  sample <- rbind(sample, paste0(str_sub(x, 1, 7), "\n", str_sub(x, 16, (nchar(sub('.*/', '', fileList[i]))-4))))
}

id <- str_sub(fileList[1], 1, 7)
file1 <- read.csv(paste0(folder, "/", fileList[1]))
file2 <- read.csv(paste0(folder, "/", fileList[2]))
file3 <- read.csv(paste0(folder, "/", fileList[3]))
file4 <- read.csv(paste0(folder, "/", fileList[4]))
file5 <- read.csv(paste0(folder, "/", fileList[5]))
file6 <- read.csv(paste0(folder, "/", fileList[6]))
file7 <- read.csv(paste0(folder, "/", fileList[7]))
file8 <- read.csv(paste0(folder, "/", fileList[8]))

GateNames <- read_csv("Panel Stainings_Tumor/_Export/GateNames_Immune2.csv")

samples <- cbind(GateNames, file1[10], file2[10])#, file3[10], file4[10], file5[10], file6[10], file7[10], file8[10])
names(samples) <- c("Population", sample[1:2])
samples[colnames(samples[2])] <- samples[,2]/samples[5,2]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[3])] <- samples[,3]/samples[5,3]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[4])] <- samples[,4]/samples[2,4]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[5])] <- samples[,5]/samples[2,5]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[6])] <- samples[,6]/samples[2,6]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[7])] <- samples[,7]/samples[2,7]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[8])] <- samples[,8]/samples[2,8]*100   # % cells of 5 immune cells / 3 live cells
samples[colnames(samples[9])] <- samples[,9]/samples[2,9]*100   # % cells of 5 immune cells / 3 live cells
# last two/six rows overwrote the samples' counts with their percentages

# 6,9,15,25,30,35,     6,9,15,25,30,35
SumPop1 <- samples[c(6,10,11,15,20,24,25,30,35),]  # immune major cell classes
#colnames(SumPop1) <- c("Population", "Baseline", "TILs")
SumPop <- as.data.table(SumPop1)
SumPop <- melt(SumPop, id.vars = c("Population"), variable.name = "Sample", value.name = "PercOfImmune")

#View(SumPop)
# cols next to each other
ggplot(data = SumPop, aes(x = Population, y = PercOfImmune, fill = Sample)) +
  geom_col(position = position_dodge(width = 0.85),
           color = "black", width = 0.65) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #legend.position = "none",
        #legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, 'cm'),
        text = element_text(size = 12)) +
  labs(title = "", y = "Cells [% of\nlive immune cells]", x = "Cell types")

ggsave(plot = last_plot(), path = folder, paste0(id, "_MajorImmuneClasses", "_TDK-vs-supernatant", ".png", sep = ""),
       width = 16, height = 8.5, units = "cm",
       bg = "transparent")

#paste("Cells in tumors")

# bars with stacked populations
ggplot(SumPop, aes(fill = Population, y = PercOfImmune, x = Sample)) + 
  geom_bar(position = "stack", stat = "identity", color = "black", width = 0.8) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        text = element_text(size = 12)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(title = paste("Cells in", "tissues"), y = "Percent of total immune cells", x = "Condition")

ggsave(plot = last_plot(), path = folder, paste0("major-immune-cells", "_VariousIDs", ".png", sep = ""),
       width = 20, height = 10, units = "cm",
       bg = "transparent")


##





############### Leftovers ######################################################################

GateNames[4,] = "Fibroblasts or\nEndothelial cells"
GateNames[8,] = "Classical\nmyeloid cells"
GateNames[13,] = "Non classical\nmyeloid cells"
GateNames[16,] = "Cytotoxic\nNK cells"
GateNames[17,] = "Pro-inflammatory\nNK cells"
GateNames[21,] = "Granulocytic\nMDSCs"
GateNames[22,] = "Monocytic\nMDSCs"
GateNames[23,] = "Neutrophils or\nMacrophages"
GateNames[30,] = "Cytotoxic\nT lymphocytes"
GateNames[31,] = "Regulatory\ncytotoxic\nT cells"
GateNames[33,] = "Tc cytotoxic\nor exhausted"
GateNames[36,] = "Conventional\nT helper cells"
GateNames[37,] = "T conv ctotoxic\nor exhausted"
GateNames[38,] = "Inflammatory\nconventional\nT cells"
GateNames[39,] = "Regulatory T\nhelper cells"
#View(GateNames)
write_excel_csv(GateNames,
                file = paste("GateNames_Immune2", ".csv", sep = ""))



Sumpop4 <- SumPop3
SumPop4 <- arrange(transform(SumPop3,
                             Population = factor(Population,levels=c("Macrophages", "Neutrophils", "Granulocytic MDSCs",
                                                                     "NKT cells", "Pro-inflammatory NK cells", "Cytotoxic NK cells",
                                                                     "cDCs", "pDCs", "Regulatory cytotoxic T cells", "RConventional T helper cells",
                                                                     "Regulatory T helper cells", "Cytotoxic T cells"))))


# ------------- T cell phenotyping panel -------------------------------------------------------

gs_pop_get_count_fast(gs[stained])$Count
#view(gs_pop_get_count_fast(gs[stained])$Count)

percentages <- data.frame(Population = rep(gs_pop_get_count_fast(gs[stained])$Population),
                          Count = rep(gs_pop_get_count_fast(gs[stained])$Count),
                          Parent = rep(gs_pop_get_count_fast(gs[stained])$ParentCount),
                          Percent = (rep(gs_pop_get_count_fast(gs[stained])$Count))/(rep(gs_pop_get_count_fast(gs[2])$ParentCount))*100)

head(percentages)
#View(percentages)

pop <- gs_pop_get_count_fast(gs[stained])$Population
pop1 <- str_extract_all(pop, "[^/](.*?)+[^/]+", simplify = TRUE)

pop1[pop1 == ""] <- NA

View(pop1)


lastValue <- function(pop1)   tail(pop1[!is.na(pop1)], 1)
#beforeLastValue <- function(pop1)  nth(pop1[!is.na(pop1)], -2L)

newGates <- apply(pop1, 1, lastValue)

newGates[19] = "Single Tumor cells"

View(newGates)

percentages["Population"] <- newGates
percentages["PercTotal"] <- percentages$Count/percentages[3,2]*100   # % cells of live cells

plot <- percentages[4:20,]



SumPop2 %>% 
  group_by(Entity) %>%
  summarise(no_samples = (length(Entity)/21))

SumPop2 %>% 
  group_by(Population) %>%
  summarise(no_samples = length(Population))



# for the Gate Names, read just the first .csv file from the folder of interest
for (i in 1:length(file_list)){
  assign(file_list[1], 
         GateNames <- read.csv(paste(folder, file_list[1], sep=''))
  )}

GateNames <- GateNames %>% rename(X1 = ?..1) %>% select(c(1:9))
test <- function(GateNames)   tail(GateNames[!is.na(GateNames)], 1)
GateNames <- apply(GateNames, 1, test)

GateNames[19] = "Other 2"
GateNames[26] = "Other 3"
GateNames[33] = "Tc ctotoxic/exhausted"
GateNames[34] = "Tc inflammatory"
GateNames[37] = "T conv ctotoxic/exhausted"
GateNames[38] = "T conv inflammatory"

GateNames <- as.data.frame(GateNames)
write_excel_csv( GateNames, file = paste("_Export/", "GateNames_Immune",
                                         ".csv", sep = ""))



# for the Gate Names, read just the first .csv file from the folder of interest
for (i in 1:length(file_list)){
  assign(file_list[1], 
         GateNames <- read.csv(paste(folder, file_list[1], sep=''))
  )}

GateNames <- GateNames %>% rename_at(1, ~"X1") %>% select(c(1:3))
test <- function(GateNames)   tail(GateNames[!is.na(GateNames)], 1)
GateNames <- apply(GateNames, 1, test)

GateNames <- as.data.frame(GateNames)
write_excel_csv( GateNames, file = paste("GateNames_Sytox",
                                         ".csv", sep = ""))







map_df_read.csv <- function(path = "_Export/Immune/", pattern = ".csv") {
  sampleFile <- read.csv(., stringsAsFactors = FALSE, sep = ",")
  cbind(GateNames, sample)
  read.csv(., stringsAsFactors = FALSE, sep = ",")[,13]
}


map_df_read.csv <- function(path = "_Export/Immune/", pattern = ".csv") {
  sampleFile <- list.files(path, pattern, full.names = TRUE) %>%
    map_df(~read.csv(.,stringsAsFactors = FALSE, sep = ","))
  cbind(GateNames, sampleFile[,13])
  read.csv(., stringsAsFactors = FALSE, sep = ",")[,13]
  }

   
files <- map_df_read.csv()
View(files)

new <- read.csv(file = "_Export/Immune/CRC_006_Immune 27-06-2021.csv")
bind_cols(GateNames, new[,-c(1:9)])

write_excel_csv(new, file = "_Export/Immune/CRC_006_Immune 27-06-2021.csv")











