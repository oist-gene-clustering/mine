##########################################################################################
#                                                                                        #
#               IMPLEMENTATION OF MAX HEAP (PRIORITY QUEUE) IN R                         #
#                                                                                        #
##########################################################################################

#########################################
# Initialization of the heap and utils. #
#########################################

#' Create heap
#'
#' Creates an empty heap.
#'
#' @return heap
#' 
#' @export
createHeap <- function() 
  return(NULL)

#' Test empty heap
#'
#' Tests if the input heap is empty.
#'
#' @param heap heap generated with this implementation
#' @return boolean
#' 
#' @export
isEmptyHeap <- function(heap) 
  return(length(heap) == 0)

#' Test if heap contains one single element
#'
#' Tests if the heap contains only one element.
#'
#' @param heap heap generated with this implementation
#' @return boolean
#' 
#' @export
oneElement <- function(heap) 
  return(length(heap) == 1)

#' Get element = key, value pair in heap
#'
#' Gets an element of the heap at a given position.
#'
#' @param pos integer giving the position of the element
#' @param heap heap generated with this implementation
#' @return element
#' 
#' @export
getElement <- function(pos, heap) 
  return(heap[[pos]])

#' Set element in heap
#'
#' Sets a element at a given position in the heap.
#' (It DOES NOT test whether the location at this position is empty)
#'
#' @param x element
#' @param pos integer giving the position to modify
#' @param heap heap generated with this implementation
#' @return modified heap
#' 
#' @export
assignValue <- function(x, pos, heap) {
  heap[[pos]] <- x
  return(heap)
}

#' Get value in heap
#'
#' Gets a value of the heap at a given position.
#'
#' @param pos integer giving the position of the element containing the value
#' @param heap heap generated with this implementation
#' @return value
#' 
#' @export
getValue <- function(pos, heap) 
  return(getElement(pos, heap)[[1]])

#' Get key in heap
#'
#' Gets a key of the heap at a given position.
#'
#' @param pos integer giving the position of the element containing the key
#' @param heap heap generated with this implementation
#' @return key
#' 
#' @export
getKey <- function(pos, heap) 
  return(getElement(pos, heap)[[2]])

#' Get root element in heap
#'
#' Gets root of the heap.
#'
#' @param heap heap generated with this implementation
#' @return root
#' 
#' @export
getMaxElement <- function(heap) 
  return(getElement(1, heap))

#' Exchange positions of elements in heap
#'
#' Exchanges the positions of two elements in a heap.
#'
#' @param father position of the first element
#' @param pos position of the second element
#' @param heap heap generated with this implementation
#' @return modified heap
#' 
#' @export
exchangeInHeap <- function(father, pos, heap) {
  mem <- getElement(father, heap)
  heap <- assignValue(getElement(pos, heap), father, heap)
  heap <- assignValue(mem, pos, heap)
  return(heap)
}

#' Add new element to heap
#'
#' Appends a new element to the heap.
#'
#' @param x element
#' @param heap heap generated with this implementation
#' @return NON CORRECT modified heap
#' 
#' @export
endQueue <- function(x, heap) 
  return(append(heap, list(x)))

#' Delete end-element of heap
#'
#' Delete the last element in the heap.
#'
#' @param heap heap generated with this implementation
#' @return modified heap
#' 
#' @export
deleteEndHeap <- function(heap) 
  return(heap[1:(length(heap)-1)])

#' Compare element keys in heap
#'
#' Compares the keys of two elements in the heap.
#'
#' @param i position of first element in heap
#' @param j position of second element in heap
#' @param heap heap generated with this implementation
#' @return 1 if key_i < key_j, -1 otherwise
#' 
#' @export
compareKey <- function(i, j, heap) 
  if (getKey(i, heap) < getKey(j, heap)) return(1) else return(-1)

##################
# Priority queue #
##################

#' Initalize queue
#'
#' Initializes the queue with a single element.
#'
#' @param x element
#' @param heap heap generated with this implementation
#' @return modified heap
#' 
#' @export
initQueue <- function(x) 
  return(list(x))

#' Push element in queue
#'
#' Pushes an element in the priority queue.
#'
#' @param value position of first element in queue
#' @param key position of second element in queue
#' @param heap heap generated with this implementation
#' @return 1 if key_i < key_j, -1 otherwise
#' 
#' @export
pushHeap <- function(value, key, heap) {
  if (isEmptyHeap(heap)) return(initQueue(list(value, key)))
  heap <- endQueue(list(value, key), heap)
  pos <- length(heap)
  father <- floor(pos/2)
  while (father > 0) {
    if (compareKey(father, pos, heap) == -1) break
    heap <- exchangeInHeap(father, pos, heap)
    pos <- father
    father <- floor(pos/2)
  }
  return(heap)
}

#' Extract max root from queue
#'
#' Extracts the root -the maximum key element- from the queue.
#'
#' @param heap heap generated with this implementation
#' @return a list containing pair root, CORRECT modified heap
#' 
#' @export
extractMaxHeap <- function(heap) {
  if (isEmptyHeap(heap)) return(NULL)
  maxi <- getMaxElement(heap)
  if (oneElement(heap)) return(list(maxi, NULL))
  pos <- 1
  heap <- exchangeInHeap(length(heap), pos, heap)
  heap <- deleteEndHeap(heap)
  child1 <- 2*pos
  child2 <- 2*pos+1
  n <- length(heap)
  while (child1 < n+1 & child2 < n+1) {
    cond1 <- (compareKey(child1, pos, heap) == 1)
    cond2 <- (compareKey(child2, pos, heap) == 1)
    if (cond1 & cond2) break
    cond3 <- (compareKey(child1, child2, heap) == 1)
    if (!cond1 & !cond3) {
      heap <- exchangeInHeap(child1, pos, heap)
      pos <- child1
      child1 <- 2*pos
      child2 <- 2*pos+1
    }
    else {
      heap <- exchangeInHeap(child2, pos, heap)
      pos <- child2
      child1 <- 2*pos
      child2 <- 2*pos+1
    }
  }
  return(list(maxi, heap))
}