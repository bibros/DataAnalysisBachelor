require(ggpubr)
require(ggsci)
require(xlsx)
library(naniar)
library(data.table)

library(openCyto)
library(ggcyto)
library(CytoML)
library(flowWorkspace)

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

# set working directory
sampleFolder <- selectDirectory(
  caption = "Select folder for working directory", label = "Select",
  path = getwd()
)
setwd(sampleFolder)

############### Get data ######################################################################

cultureSuccess <- read_excel(path = "ALI_organoid_CultureSuccess.xlsx", sheet = "Culture time success")
cultureContami <- read_excel(path = "ALI_organoid_CultureSuccess.xlsx", sheet = "Culture time contaminated")
cultureNoGrowt <- read_excel(path = "ALI_organoid_CultureSuccess.xlsx", sheet = "Culture time not grown")

list.df <- list(data.frame(cultureSuccess[,2]), data.frame(cultureContami[,2]), data.frame(cultureNoGrowt[,2]))
max.rows <- max(unlist(lapply(list.df, nrow), use.names = F))

list.df <- lapply(list.df, function(x) {
  na.count <- max.rows - nrow(x)
  if (na.count > 0L) {
    na.dm <- matrix(NA, na.count, ncol(x))
    colnames(na.dm) <- colnames(x)
    rbind(x, na.dm)
  } else {
    x
  }
})

cultureTime <- do.call(cbind, list.df)
colnames(cultureTime) <- c("Successful culture", "Contaminated", "Did not grow")
cultureTime2 <- cultureTime

for (i in 1:3){
  for(j in 1:length(cultureTime[,i])){
    if(!is.na(cultureTime[j,i])){
      cultureTime2[j,i] <- sub("\\_.*", "", cultureTime[j,i])
    }
  }
}


de <- merge(table(cultureTime2$`Successful culture`),
            table(cultureTime2$Contaminated), by = 1, all=TRUE)
fr <- merge(de,
            table(cultureTime2$`Did not grow`), by = 1, all=TRUE)
colnames(fr) <- c("Entity", "Success", "Contaminated", "No growth")

#fr %>% arrange(Entity)
x <- c("BLD", "CRC", "GBL", "HNO", "LUN", "LVR", "MAM", "MEL", "MNG",
       "OES", "OVA", "PAN", "RCC", "SCL", "SMI", "STM", "TES", "THY")

fr <- fr %>% slice(match(x, Entity))
fr <- as.data.table(fr)
fr$Entity <- factor(fr$Entity, levels = fr$Entity)

data <- melt(fr, id.vars = "Entity", variable.name = "Growth condition", value.name = "Sample numbers")

# Stacked
ggplot(data, aes(fill = Entity, y = `Sample numbers`, x = `Growth condition`)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values = c("#F33C6C", "#EA830E", "#D9A628", "#BED228", "#9CA739", "#6EB043", "#009A4B", 
                               "#2BB45F", "#2FA18B", "#28BCB1", "#1AABC3", "#26C3E9", "#4A9FE2", "#837CC8", 
                               "#A283BC", "#B978B3", "#D16FAC", "#EA68A7", "#F0B097")) +
  guides(fill = guide_legend(nrow = 3)) +
  labs(title = "", y = "Number of tumor samples", x = "Condition") +
  coord_flip()

ggsave(plot = last_plot(), paste("Tumor_culture_success2", ".png", sep = ""),
       width = 18, height = 8, units = "cm",
       bg = "transparent")



# scale_fill_manual(values = c("#FF0000", "#FF6200", "#FFAA00", "#FFFF00", "#99FF00", "#00E500", "#009600", 
                              #"#006C00", "#006C7B", "#0000FF", "#006CE9", "#00A6FF", "#9AA4FF", "#FFA4FF", 
                              #"#FF59FF", "#B800FF", "#7100FF", "#93008E", "#930016")) +