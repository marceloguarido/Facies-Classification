#' Function test
plotwelllogs <- function(datap, wname, dz, wlogs, facies, fcolor, fnames){
  # Used libraries
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(ggpubr)
  
  # Selecting only the data for the chosen well
  datat <- datap[datap$Well.Name == wname,]
  datat[ , !(names(datat) %in% "Well.Name")]
  
  # Creating a column of ones to use on the facies plotting
  datat$temp = rep(1,nrow(datat))
  
  # Organizing the data for the plot
  df1 <- melt(datat[,c(dz,wlogs)], id = c(dz))
  df2 <- datat[,c(dz,"temp",facies)]
  colnames(df1) <- c("Depth", "variable", "Measurements")
  colnames(df2) <- c("Depth", "temp", "Facies")
  
  # Using the melted df for the facet_grid plot (subploting all the well logs)
  p1 <- ggplot(df1, aes(Depth, Measurements, group=variable, color=variable)) + 
    facet_grid(. ~ variable, scales = "free_x") +
    geom_line() +
    coord_flip() +
    scale_x_reverse() +
    theme(legend.position = "none",
          strip.background =element_rect(fill = "gray55"),
          strip.text = element_text(face = "bold", color = "white", size = 12),
          axis.text.x = element_text(angle = 45, size = 12, face = "bold", hjust = 1),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)
          )
  
  # Creating the tile plot for the facies classification plot
  p2 <- ggplot(df2, aes(y = Depth, x = temp)) + 
    geom_tile(aes(fill = factor(Facies))) +
    scale_fill_manual(values = facies_colors,
                      labels = fnames,
                      name = "Facies") +
    labs(title = "Facies") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(color = "white"),
          axis.text.x = element_text(color = "white", size = 12, angle = 45, face = "bold"),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(color = "white", face = "bold")
    ) +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 1.65))
  
  # Using the well name as the plot title
  tname = text_grob(wname,
                    face = "bold.italic", 
                    color = "steelblue",
                    size = 22)
  # Legend name
  lname = text_grob("Facies",
                    face = "bold",
                    color = "darkblue",
                    size = 16,
                    rot = 270)
  
  # Griding the facet_plot and the tile plot in one single figure
  grid.arrange(p1,p2,
               nrow = 1,
               ncol = 2,
               widths = c(length(wlogs),1.2),
               top = tname,
               right = lname)
}


##### Function to plot the Confusion Matrix in a more fashion way
plotConfusionMatrix <- function(CMtable, relScale = TRUE){
  # Used libraries
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(ggpubr)
  library(caret)
  
  if(relScale == TRUE){
    xtab = round(CMtable/rowSums(CMtable)*1000)/1000
    limits = c(0,1)
    tname = c("Confusion Matrix: Relative Values")
  }
  else {
    xtab = CMtable
    limits = c(0,max(CMtable))
    tname = c("Confusion Matrix: Absolute Values")
  }
  
  meltedCM <- melt(xtab)
  colnames(meltedCM) <- c("Predictions", "True", "value")

  ggplot(data = meltedCM, aes(x = Predictions, y = True, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "lightblue", high = "orchid", mid = "snow",
                         midpoint = 0,
                         limit = limits, space = "Lab") +
    scale_x_discrete(limits = rev(levels(meltedCM$Predictions))) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                     size = 8, hjust = 0.5),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 11)) +
    geom_text(aes(Predictions, True, label = value), color = "black", size = 3) +
    ggtitle(tname) +
    coord_flip() +
    scale_y_discrete(position = "right") 
    # coord_fixed()
}


##### Function to convert cartesian coordinates to polar coordinates
cart2polwells <- function(data_in, features_wells){
  data_polar = data_in
  fea_red = features_wells
  for (fea1 in features_wells){
    fea_red <- fea_red[-1]
    for (fea2 in fea_red){
      # x = data_in[,fea1]
      # y = data_in[,fea2]
      x = data_in[,fea1] - mean(data_in[,fea1])
      y = data_in[,fea2] - mean(data_in[,fea2])
      # x = x / max(abs(x))
      # y = y / max(abs(y))
      data_polar[,c(paste(fea1,"_",fea2,"_rho", sep = ""))] = sqrt(x^2 + y^2)
      data_polar[,c(paste(fea1,"_",fea2,"_phi", sep = ""))] = atan(y/x)
    }
  }
  return(data_polar)
}



### Function to return the legend for plots
getLegend <- function(data, wname, dz, wlogs, facies, facies_colors, fnames){
datat <- data[data$Well.Name == wname,]
datat[ , !(names(data) %in% "Well.Name")]

# Creating a column of ones to use on the facies plotting
datat$temp = rep(1,nrow(datat))

# Organizing the data for the plot
df1 <- melt(datat[,c(dz,wlogs)], id = c(dz))
df2 <- datat[,c(dz,"temp",facies)]
colnames(df1) <- c("Depth", "variable", "Measurements")
colnames(df2) <- c("Depth", "temp", "Facies")
p2 <- ggplot(df2, aes(y = Depth, x = temp)) + 
  geom_tile(aes(fill = factor(Facies))) +
  scale_fill_manual(values = facies_colors,
                    labels = fnames,
                    name = "Facies") +
  labs(title = "Facies") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white", size = 6.5, angle = 45, face = "bold"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_text(size = 10, face = "bold", color = "darkblue"),
        legend.text = element_text(size = 7),
        plot.title = element_text(color = "white", face = "bold")
  ) +
  guides(fill = guide_legend(keywidth = 0.7, keyheight = 1.65))

pright <- grab_legend(p2)
return(pright)
}

#### Compute the gradient over depth of the features
features_gradient <- function(df, features, wnames, depth = "Depth", wncol = "Well.Name"){
  for (wname in wnames){
    dfT <- df[df[ ,wncol] == wname, ]
    
    for (feature in features){
      df[which(df[ ,wncol] %in% c(wname)), paste0(feature, "_grad")] = c(0, diff(as.matrix(dfT[ ,feature]))/diff(as.matrix(dfT[ , depth])))
    }
    
  }
  
  return(df)
}

#### Function to create confusion matrix with no need of the factors be of the same size
createConfusionMatrix <- function(act, pred) {
  numClasses <- max(act, pred)
  pred <- pred[order(act)]
  act  <- act[order(act)]
  sapply(split(pred, act), tabulate, nbins=numClasses)
}

### Function to plot the well logs with the predictions
plotwellPred <- function(datap, wname, dz, wlogs, facies, classy, fcolor, fnames){
  # Used libraries
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
  library(ggpubr)
  
  # Selecting only the data for the chosen well
  datat <- datap[datap$Well.Name == wname,]
  datat[ , !(names(datat) %in% "Well.Name")]
  
  # Creating a column of ones to use on the facies plotting
  datat$temp = rep(1,nrow(datat))
  
  # Organizing the data for the plot
  df1 <- melt(datat[,c(dz,wlogs)], id = c(dz))
  df2 <- datat[,c(dz,"temp",facies)]
  df3 <- datat[,c(dz,"temp",classy)]
  colnames(df1) <- c("Depth", "variable", "Measurements")
  colnames(df2) <- c("Depth", "True", "Facies")
  colnames(df3) <- c("Depth", "Pred", "Predictions")
  
  # Using the melted df for the facet_grid plot (subploting all the well logs)
  p1 <- ggplot(df1, aes(Depth, Measurements, group=variable, color=variable)) + 
    facet_grid(. ~ variable, scales = "free_x") +
    geom_line() +
    coord_flip() +
    scale_x_reverse() +
    theme(legend.position = "none",
          strip.background =element_rect(fill = "gray55"),
          strip.text = element_text(face = "bold", color = "white", size = 12),
          axis.text.x = element_text(angle = 45, size = 12, face = "bold", hjust = 1),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)
    )
  
  # Creating the tile plot for the facies classification true plot
  p2 <- ggplot(df2, aes(y = Depth, x = True)) + 
    geom_tile(aes(fill = factor(Facies)), show.legend = F) +
    scale_fill_manual(values = facies_colors,
                      labels = element_blank(),
                      name = "Facies") +
    labs(title = "True") +
    ggtitle("True") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(color = "white"),
          axis.text.x = element_text(color = "white", size = 12, angle = 45, face = "bold"),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(color = "steelblue", face = "bold", hjust = 0.5, size = 16)
    ) +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 1.65))
  
  # Creating the tile plot for the facies classification predictions plot
  p3 <- ggplot(df3, aes(y = Depth, x = Pred)) + 
    geom_tile(aes(fill = factor(Predictions))) +
    scale_fill_manual(values = facies_colors,
                      labels = fnames,
                      name = "Facies") +
    labs(title = "Pred") +
    ggtitle("Pred") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(color = "white"),
          axis.text.x = element_text(color = "white", size = 12, angle = 45, face = "bold"),
          axis.ticks.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(color = "steelblue", face = "bold", hjust = 0.5, size = 16)
    ) +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 1.65))
  
  # Using the well name as the plot title
  tname = text_grob(wname,
                    face = "bold.italic", 
                    color = "steelblue",
                    size = 22)
  # Legend name
  lname = text_grob("Facies",
                    face = "bold",
                    color = "darkblue",
                    size = 16,
                    rot = 270)
  
  # Griding the facet_plot and the tile plot in one single figure
  grid.arrange(p1,p2,p3,
               nrow = 1,
               ncol = 3,
               widths = c(length(wlogs),.7,1.2),
               top = tname,
               right = lname)
}
