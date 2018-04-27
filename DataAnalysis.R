getwd()
library(ggplot2)


file1 <- read.csv("/home/freek/Introgression/Data/25-04-2018-14-47-56/Data.csv", header = T)
file2 <- read.csv("/home/freek/Introgression/Data/25-04-2018-14-53-15/Data.csv", header = T)
file3 <- read.csv("/home/freek/Introgression/Data/25-04-2018-14-58-36/Data.csv", header = T)
file4 <- read.csv("/home/freek/Introgression/Data/25-04-2018-15-15-49/Data.csv", header = T)
file5 <- read.csv("/home/freek/Introgression/Data/25-04-2018-15-21-18/Data.csv", header = T)

file1$rescue <- file1$AVG_p/file1$AVG_q
file2$rescue <- file2$AVG_p/file2$AVG_q
file3$rescue <- file3$AVG_p/file3$AVG_q
file4$rescue <- file4$AVG_p/file4$AVG_q
file5$rescue <- file5$AVG_p/file5$AVG_q

file1$pnq <- file1$AVG_p+file1$AVG_q
file2$pnq <- file2$AVG_p+file2$AVG_q
file3$pnq <- file3$AVG_p+file3$AVG_q
file4$pnq <- file4$AVG_p+file4$AVG_q
file5$pnq <- file5$AVG_p+file5$AVG_q

ggplot(file1, aes(x=Generation, y =rescue)) +
    geom_line(aes(x=Generation, y =rescue), color = 'red') + 
    geom_ribbon(aes(ymin=rescue-VAR_p,ymax=rescue+VAR_p), color = 'red') + 
    geom_line(data = file2, aes(x=Generation, y =rescue), color = 'blue') + 
    geom_ribbon(data = file2, aes(ymin=rescue-VAR_p,ymax=rescue+VAR_p), color = 'blue') +
    geom_line(data = file3, aes(x=Generation, y =rescue), color = 'green') + 
    geom_ribbon(data = file3, aes(ymin=rescue-VAR_p,ymax=rescue+VAR_p), color = 'green') +
    geom_line(data = file4, aes(x=Generation, y =rescue), color = 'purple') + 
    geom_ribbon(data = file4, aes(ymin=rescue-VAR_p,ymax=rescue+VAR_p), color = 'purple') +
    geom_line(data = file5, aes(x=Generation, y =rescue)) + 
    geom_ribbon(data = file5, aes(ymin=rescue-VAR_p,ymax=rescue+VAR_p))


plot(file1$rescue~file1$Generation, ylim = c(0,0.05), type = 'l', col = 'red', lwd=5, ylab = "p / q", xlab = "Generation")
par(new = T)
plot(file2$rescue~file1$Generation, ylim = c(0,0.05), type = 'l', col = 'green', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file3$rescue~file1$Generation, ylim = c(0,0.05), type = 'l', col = 'purple', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file4$rescue~file1$Generation, ylim = c(0,0.05), type = 'l', col = 'blue', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file5$rescue~file1$Generation, ylim = c(0,0.05), type = 'l', col = 'orange', lwd=5, ylab = " ", xlab = " ")

plot(file1$pnq~file1$Generation, ylim = c(0,2), type = 'l', col = 'red', lwd=5, ylab = "p+q", xlab = "Generation")
par(new = T)
plot(file2$pnq~file1$Generation, ylim = c(0,2), type = 'l', col = 'green', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file3$pnq~file1$Generation, ylim = c(0,2), type = 'l', col = 'purple', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file4$pnq~file1$Generation, ylim = c(0,2), type = 'l', col = 'blue', lwd=5, ylab = " ", xlab = " ")
par(new = T)
plot(file5$pnq~file1$Generation, ylim = c(0,2), type = 'l', col = 'orange', lwd=5, ylab = " ", xlab = " ")
    

file1 <- read.csv("/home/freek/Introgression/Data/26-04-2018-15-35-54/AlleleF_mean.csv", header = F)

filematrix <- data.matrix(file1)

plot(filematrix[50,], ylim = c(0,1))


## START NEW PROJECT:

list.data<-list()
list.param<-list()
rec<-c()
dist <- c()
nin1 <- c()
introgressed<-c()

filenames <- list.files(paste(getwd(),"/Data", sep = ''))

for(i in 1: length(filenames)){
    list.data[[i]] <- read.csv(paste(paste(paste(getwd(),"/Data",sep = ''), (filenames[i]),sep='/'),"/Data.csv", sep =''),  header = TRUE, sep = ',')
    list.param[[i]] <- read.csv(paste(paste(paste(getwd(),"/Data",sep = ''), (filenames[i]),sep='/'),"/Parameters.csv", sep =''), header = FALSE) 
    rec[i] <- list.param[[i]]$V2[11]
    nin1[i] <- list.param[[i]]$V2[5]
    dist[i] <- list.param[[i]]$V2[14]
    introgressed[i] <- list.data[[i]]$Generation[length(list.data[[i]]$Introgressed)]
}

data <- data.frame(rec,nin1,dist,introgressed)

ggplot(data, aes(x=rec,y=introgressed, color = dist))+ geom_point()
