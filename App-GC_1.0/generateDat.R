generateDat <- function(n)
{
  dat = data.frame()
  for(i in 1:n)
  dat = rbind(dat,mtcars[i,])
  dat
}

generateDat2 <- function(coordinate)
{
  dat = data.frame()
  for(i in 1:4)
    dat = rbind(dat,mtcars[i,])
  dat = cbind(dat,unlist(coordinate))
  dat
}