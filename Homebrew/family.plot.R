#pedigree.plot

#inputs a COLONY file with a column for offspring, moms, and dads
#Does NOT include a cohort column
#should accept a BestCluster or BestColony file

family.plot <- function(family,title = "Pedigree Plot"){
  family <- family[order(family$MotherID,family$FatherID),]
  moms <- unique(family$MotherID)
  dads <- unique(family$FatherID)

  momdots <- data.frame(ID = unique(family$MotherID), y = seq(from=1,to=length(family$OffspringID),length.out=length(moms)))
  daddots <- data.frame(ID = unique(family$FatherID), y = seq(from=1,to=length(family$OffspringID),length.out=length(dads)))
  
  #make the plot points
  plot(x = c(1,2,3), 
       y = c(0, length(family$OffspringID), length(family$OffspringID)), 
       pch = "", axes = FALSE, 
       xlab = "", ylab = "",
       main = title)
  points(x = rep(2,length(family$OffspringID)), y = c(1:length(family$OffspringID)), pch = "-", cex = 0.25)
  points(x = rep(1,length(moms)), y = momdots$y, pch = 19, cex = 0.5)
  points(x = rep(3,length(dads)), y = daddots$y, pch = 19, cex = 0.5)
  
  #make the lines in the plot
  for(r in 1:length(family$OffspringID)) {
    #mom line
    lines(x=c(2,1), y = c(r, momdots[which(momdots$ID == family$MotherID[r]),]$y), lwd = 0.5,col = "grey")
    
    #dad line
    lines(x =c(2,3), y = c(r, daddots[which(daddots$ID == family$FatherID[r]),]$y), lwd = 0.5,col = "grey")
  }
  
  #adding labels
  mtext("Parent 1", side = 3, line = 0, at = 1, cex = 0.7)
  mtext("Parent 2", side = 3, line = 0, at = 3, cex = 0.7)
  mtext("Offspring", side = 3, line = 0, at = 2, cex = 0.7)
}
