# Functions for plotting cells in four field Neftel space

# Calculating D score
calculate_D <- function(row) {
    D <- max(as.numeric(row[c("OPC_State","NPC1_state","NPC2_state")]))-max(as.numeric(row[c("AC_state","MES1_state","MES2_state")]))
    return(D)
}

# OPC/NPC cells have (D>0) and x axis value is defined using following formula: log2(|SCopc – SCnpc|+1)
# AC/MES cells have (D<0) and x axis value is defined using following formula: log2(|SCac – SCmes|+1).
# AC and OPC x axis values are multiplied with -1 to separate the to their own quarter of the plot.

calculate_x<-function(row){
  if(as.numeric(row["D"])>0){ # cell is OPC/NPC
    x<-log2(abs(as.numeric(row["OPC_state"]) - max(as.numeric(row[c("NPC1_state","NPC2_state")])))+1)
  }else{ # cell is AC/MES
    x<-log2(abs(as.numeric(row["AC_state"]) - max(as.numeric(row[c("MES1_state","MES2_state")])))+1)
  }
  if(row[["neftel_state"]] %in% c("OPC_state","AC_state")){
    print(paste("Flipping sign for:", row[["neftel_state"]])) 
    x<-x*-1
  }
  return(x)
}

# Four field scatter plot of the cells
plot_neftel_space = function(metadata, column) {
    if (is.character(metadata[[column]])) {
        metadata[[column]] = factor(metadata[[column]], levels=unique(metadata[[column]]))
    } 
     quad_labels <- data.frame(x = c(0.65, -0.65, -0.65, 0.65), y = c(-0.65, -0.65, 0.65, 0.65),
                              label = c("MES-like", "AC-like", "OPC-like", "NPC-like"))

    p = ggplot(metadata, aes(x=x, y=D, color=.data[[column]], alpha=.data[[column]])) + 
        geom_point(size=1, shape=16)+
        geom_text(data=quad_labels, aes(x=x, y=y, label=label), inherit.aes=F, color="black") +
        xlab("log2(|SC1 - SC2| + 1)")+
        ylab("log2(|SC1 - SC2| + 1)")+
        theme_classic()+
        geom_vline(xintercept = 0,linetype="dashed",col="grey")+
        geom_hline(yintercept = 0,linetype="dashed",col="grey")+
        xlim(c(-0.7,0.7))+ylim(c(-0.7,0.7))
    return(p)
}

