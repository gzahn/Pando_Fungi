find_gps_dists <- 
function(points1,points2){
  
  # tests
  stopifnot(class(points1) == "data.frame")
  stopifnot(class(points2) == "data.frame")
  
  if(ncol(points1) != 2 | ncol(points2) != 2){
    stop("data frames must have 2 columns only")
  }
  
  if(!apply(points1,2,class) %>% unique() %in% c("numeric","integer")){
    stop("columns must be numeric; col1=longitude,col2=latitude")
  }
  
  if(!apply(points2,2,class) %>% unique() %in% c("numeric","integer")){
    stop("columns must be numeric; col1=longitude,col2=latitude")
  }
  
  # actual function
  mylist <- list()
  for(i in 1:nrow(points1)){
    mylist[[i]] <- points1[i,] %>% unlist
  }
  
  distances <- list()
  for(i in 1:nrow(mypoints)){
    
    mydistfunction <- function(x){geosphere::distHaversine(x,points2[i,])}
    colname <- paste0("dist_to_",i)
    distances[[colname]] <- map_dbl(mylist, mydistfunction)
  }
  return(as.data.frame(distances))
}
