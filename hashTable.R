##############################################################
#                                                            #
#               HASH TABLE IMPLEMENTATION                    #
#                                                            #
##############################################################

#########################################
# Initialization of Hash Table & utils. #
#########################################
## Most of the functions are actually wrappers for R environments                    ##
# https://blog.dominodatalab.com/a-quick-benchmark-of-hashtable-implementations-in-r/ #

#' Create Hash Table
#'
#' Creates a Hash Table as a R environment.
#'
#' @return hash table
#' 
#' @export
createHT <- function() return(new.env(hash=T))

#' Create key for hash table
#'
#' Creates the key for the hash table from the key provided.
#'
#' @param key vector of two integers
#' @return real hash table key
#' 
#' @export
getName <- function(key) return(paste0("hashTableArray[", paste(key[1],key[2],sep=","), "]"))

#' Compute value and store it in the Hash Table
#'
#' Computes a value and stores it in the Hash Table.
#'
#' @param key key associated with value
#' @param ht hash table
#' @param fn function with a single argument @key
#' @return value
#' 
#' @export
computeNStore <- function(key, ht, fn=NULL) {
  value <- fn(key)
  assign(getName(key), value, envir=ht)
  return(value)
}

###################
# Main operations #
###################

## Values (that should be valid regular R objects) are stored with an unique key, which is a vector of two integers ##

#' Get value from the Hash Table
#'
#' Gets value from the Hash Table.
#'
#' @param key key associated with value
#' @param ht hash table
#' @param fn function with a single argument @key
#' @return value
#' 
#' @export
getHT <- function(key, ht, fn=NULL)
  return(tryCatch(get(getName(key), envir=ht), 
           error=function(e) return(computeNStore(key, ht, fn))
   ))
  