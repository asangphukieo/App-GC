countdown <- function(from)
{
  print(from)
  while(from!=0)
  {
    Sys.sleep(1)
    from <- from - 1
    print(from)
  }
  x=c(1:3)
  x
}


