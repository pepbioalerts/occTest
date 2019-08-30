
loadSettings <- function (x){
  for(i in 1:length(x)) {

    a = x
    
    for (ii in 1:length(a)){ 
      assign(names(a)[ii], a[[ii]],inherits = T,envir = parent.env())
    }
    
    
    
  }
}




  
