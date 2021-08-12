setwd("~/Uni/6) Bachelorarbeit/6 Data")

install.packages("hrbrthemes")
install.packages("viridis")
install.packages("packcircles")
install.packages("RColorBrewer")
install.packages("ggsci")
install.packages("pals")
install.packages("readxl")

require(tidyverse)
require(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(packcircles)
library(RColorBrewer)
library(ggsci)
require(pals)
require(readxl)
library(scales)

data1 <- read_excel("summary.xlsx") %>% select(c(Tumor, count))
head(data1)


############ Pie chart ###############################################################################

data1 <- as.data.frame(data1)

data <-
  arrange(data1, desc(Tumor)) %>%
  mutate(percent = count / sum(data1$count) *100) %>%
  mutate(ypos = (cumsum(percent) - 0.5*percent))

ggplot(data, aes(x = "", y = percent, fill = Tumor)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = ypos, label = Tumor), color = "black", size = 5)+
  theme_void()

ggplot(data, aes(x = 2, y = percent, fill = Tumor)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0)+
  geom_text(aes(y = ypos, label = Tumor), color = "black", size = 5)+
  theme_void()+
  xlim(0.75, 2.5)


############ Bar plot ################################################################################

ggplot(data1, aes(x = reorder(Tumor, -count), y = count, fill = Tumor)) +
  geom_bar(stat = "Identity", color = "black") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(size = 14),
        legend.position = "none",
        text = element_text(size = 15)) +
  labs(title = "Sample numbers of different tumor entities",
       y = "Count of tumors", x = "Tumor entities")
  
ggsave(plot = last_plot(), paste("Tumor_Entities_bar_", format(Sys.time(), "%Y-%m-%d"),
                                 ".png", sep=""), 
       width = 20, height = 10, units = "cm",
       bg = "transparent")
ggsave(plot = last_plot(), paste("Tumor_Entities_bar_", format(Sys.time(), "%Y-%m-%d"),
                                 ".pdf", sep = ""),
       width = 21, height = 29, units = "cm",
       bg = "transparent")
  




############ Bubble plot #############################################################################


# Generate the layout. This function return a dataframe with one line per bubble. 
# It gives its center (x and y) and its radius, proportional of the value
packing <- circleProgressiveLayout(data1$count, sizetype='area')

# We can add these packing information to the initial data frame
data1 <- cbind(data1, packing)
write.csv(data1,"summary_new.csv", row.names = FALSE)

# Check that radius is proportional to value. We don't want a linear relationship, since it is the AREA that must be proportionnal to the value
# plot(data$radius, data$value)

# The next step is to go from one center + a radius to the coordinates of a circle that
# is drawn by a multitude of straight lines.
dat.gg <- circleLayoutVertices(packing, npoints=500)

# Make the plot
ggplot() + 
    # Make the bubbles
  geom_polygon(data = dat.gg, aes(x=x, y=y, group = id, fill=as.factor(id)),
               color = "black", show.legend = F) +
    # Add text in the center of each bubble + control its size
  geom_text(data = data1, aes(x=x, y=y, label = Tumor, size = count)) +
  scale_size_continuous(range = c(5,9)) +
  
  # General theme:
  theme_void() +
  scale_fill_manual(values = cubicyf(nrow(data1))) +
  #scale_fill_manual(values = colorRampPalette(pal.bands(8, "stepped"))(15)) +
  theme(legend.position="none") +
  coord_equal()


ggsave(plot = last_plot(), paste("Tumor_Entities_bubble_", format(Sys.time(), "%d-%m-%Y"),
                                 ".png", sep=""), 
       width = 20, height = 20, units = "cm",
       bg = "transparent")
ggsave(plot = last_plot(), paste("Tumor_Entities_bubble_", format(Sys.time(), "%d-%m-%Y"),
                                 ".pdf", sep = ""),
       width = 21, height = 29, units = "cm",
       bg = "transparent")

#kovesi.rainbow
#stepped3
#cubicyf
#brewer.gnbu

# look here for color palettes of the pal package:
# https://cran.r-project.org/web/packages/pals/vignettes/pals_examples.html
pal.bands(coolwarm, parula, ocean.haline, brewer.blues, cubicl, kovesi.rainbow, ocean.phase, brewer.paired(12), stepped, brewer.seqseq2,
          main="Colormap suggestions")
pal.channels(parula, main = "parula")

display.brewer.all()






# my own tries with given x/y values (which were by the way the same as the function above created...)

# Most basic bubble plot
data2 <- data %>% #filter(Tumor == "RCC") %>%
  mutate(
    x = replace(x, Tumor, c(-1.2, 1.25, 1.45, 0.48, 3.27,
                     -1.18, -3.7, -4, -2.15, -2.137,
                     0.75, 2.6, 4.2, 4.472, -2.138)),
    y = replace(y, Tumor, c(-0.3, 0.17, -2.75, 2.75, 1.2,
                     -2.8, -2.2, 0.7, 2.204, 5.648,
                     5.4, 4.8, 3, 1.897, -4.339)))

#View(data2)
ggplot(data2, aes(x=x, y=y, size = count, color = Tumor, stroke = 1.5)) +
  geom_point(show.legend = F) +
  geom_text(aes(x=x, y=y, label = Tumor), colour = I(alpha("black")), size = 3) +
  xlim(-8, 9) +  ylim(-7, 8) +
 # scale_linetype_manual(100) +
  scale_size(range = c(13,35)) +
  scale_fill_viridis(discrete=TRUE, guide=FALSE, option="A") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

