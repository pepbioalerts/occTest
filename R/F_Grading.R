### FUNCTIONS FOR BUILDING QUALITY GRADING 

.lazylogic <- function (e){o <- eval(parse(text=e)); ifelse(is.na(o),F,o) }

.make.qgrading.system <- function (codes=codes.list,
                                  conditionals.codes=conditionals.codes.list, 
                                  descriptors.codes=descriptors.codes.list, verbose=F){
  
  
  #making sure everything is right
  stopifnot (length(codes) == length(conditionals.codes))
  stopifnot (length(codes) == length(descriptors.codes))
  stopifnot (is.list(conditionals.codes))
  stopifnot (is.list(descriptors.codes))
  stopifnot (is.character(codes))
  stopifnot (is.na (conditionals.codes[[length(conditionals.codes)]]) )
  stopifnot (is.na (descriptors.codes[[length(descriptors.codes)]]) )
  
  l <- as.list (rep(NA, length(conditionals.codes))) 
  names(l) <- codes
  
  for (i in 1:length(l)) { 
    
    l[[i]] <- list (Conditional = conditionals.codes[[i]][1], Descriptor = descriptors.codes[[i]] )  
    
  }
  
  
  return (l)
  
  
  
  
  
  
  
  
}


.make.qgrading <- function (incase, .grading.scheme, verbose=F){
  attach(incase,warn.conflicts = F)

  for (i in 1:length(.grading.scheme)){

    #check if we need to continue because it is the maximum quality
    if ( i==length(.grading.scheme) ) { 
      qlabel <- names (.grading.scheme)[i]  
      output.descriptor <- NA
      detach(incase)
      break
    }
    
    
    grading.params <- .grading.scheme[[i]]
    result.grading <-  eval(parse(text=grading.params$Conditional))
    
    #check if we need to continue because it did not meet the grading classification
    if (result.grading==F) {next()}
    
    
    #get quality label
    if(result.grading==T) {qlabel <- names (.grading.scheme)[i]}
    #get descriptor attributes
    if(result.grading==T){
      
      output.descriptor <- lapply(grading.params$Descriptor,function (descriptor) {
        
        o <- try (eval(parse(text=descriptor)),silent = T)
        
        if (class (o) =='try-error') {return (descriptor)} else {return(o)}
        
      })
      output.descriptor <- .subsetlist.nonNULL(output.descriptor)
      output.descriptor <- unlist (output.descriptor)
      output.descriptor <- paste(output.descriptor, collapse = '-')
    
    }
    #stop it if already found the thing
    if(result.grading==T) {detach(incase); break}
    
  }
  
  output.grading <- data.frame (quality.grade=qlabel,quality.descriptor=output.descriptor) 
  
  return(output.grading)
  
}



.qgrade.data.frame <-  function (df,grading.scheme, verbose=F){
  
  if(verbose) print('Performing QAQC')
  q.out <-  lapply (1:nrow (df), function (r){
    out <- occProfileR:::.make.qgrading (incase = df[r,] , .grading.scheme =grading.scheme )
    return (out)
  })
  
  if (verbose) print('Coercing to data.frame')
  q.out <- do.call(rbind,q.out)
  
  return (q.out)
  
}



