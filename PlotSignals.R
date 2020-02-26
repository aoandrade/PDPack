library('dygraphs')


# This utility files provides functions to plot signals from a dataframe



plotAccelerometer <- function(df){
  gg1<-qplot(x=df$X.Time., df$X.A1.X.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.X")
  gg2<-qplot(x=df$X.Time., df$X.A1.Y.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Y")
  gg3<-qplot(x=df$X.Time., df$X.A1.Z.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Z")
  gg4<-qplot(x=df$X.Time., df$X.A2.X.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.X")
  gg5<-qplot(x=df$X.Time., df$X.A2.Y.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.Y")
  gg6<-qplot(x=df$X.Time., df$X.A2.Z.)+
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.Z")
  # gg7<-qplot(x=df$X.Time., df$X.A3.X.)+
  #   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.X")
  # gg8<-qplot(x=df$X.Time., df$X.A3.Y.)+
  #   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.Y")
  # gg9<-qplot(x=df$X.Time., df$X.A3.Z.)+
  #   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.Z")
  
  
  grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = "") #, gg7, gg8, gg9, nrow = 3, ncol=3, top = "")
  
}



ggplotAccelerometer <- function(df){
  gg1<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A1.X.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.X") +
    geom_line()
  gg2<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A1.Y.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Y") +
    geom_line()
  gg3<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A1.Z.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Z") +
    geom_line()
  gg4<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A2.X.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.X") +
    geom_line()
  gg5<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A2.Y.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.Y") +
    geom_line()
  gg6<-ggplot(data = df, mapping = aes(x = df$X.Time., y = df$X.A2.Z.)) +
    xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A2.Z") +
    geom_line()
  grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = "") 
  
}
