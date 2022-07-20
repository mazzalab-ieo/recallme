#function for computing the Jaccard index
jaccard <- function(x, y) {
  intersection = length(intersect(x, y))
  union = length(x) + length(y) - intersection
  return (intersection/union)
}

#function for computing the cosine similarity
cosine_similarity <- function(x,y) {
  xtemp = sample(x, length(y) )
  return((xtemp%*%y) / ( sqrt(sum(xtemp^2)) * sqrt(sum(y^2)) )  )
}

