
#' Create value of items in a collection
#'
#' @param collection A data frame containing a collection.
#' @param dependents A list of scoring systems, where a scoring system is a list of:
#' - `$column` the column in the `collection` we are generating a score from.
#' - `$weight` the weight of this score (multiplier), relative to other columns.
#' - `$details` a data frame detailing the score for particular values found in the column. Where the values to match to the column are found in `names` and the score is contained in `value`.
#' - `$type` The method for matching/finding entries in the column. Allowable values are `'Exact'`, `'Grepl'`, `'Number range'`  and `'Number range (int)'`.
#'
#' @details
#' This function generates a value for each item in a collection given a users scoring system. The scoring system is provided via the input `dependents`. Please see the vignette `Adding value to items in a collection` for further documentation.
#'
#' @return A list of length two containing:
#' - $total_value a vector of length `nrow(collection)` containing the value of each item in the collection.
#' - $value_breakdown a data frame with `nrow(collection)` rows and `length(dependents)` columns giving the value of each dependent in corresponding columns.
#' @export
#'
value_of_items <- function(collection, dependents = NA){

  ### 1) Checks.
  if(is.na(dependents[1])){
    stop('Require dependents!')
  }
  columns = unlist(lapply(dependents,function(x){x$column}))
  if(!all(columns %in% names(collection))){
    stop('Error! Not all columns in `dependents` are contained within `collection`')
  }


  ### 2) Create value_start vector of 0s.
  value_start = rep(0, nrow(collection))

  ### 3) Loop over all dependents and create the value vector for each.
  values_each_dependent = data.frame(lapply(dependents, function(dependent){
    value = value_start

    quantity = collection[[dependent$column]]

    if(dependent$type == 'Exact'){
      # In this case we use '==' to check quantity aganist dependent$details$names
      for(i in 1:nrow(dependent$details)){
        value[quantity == dependent$details$names[i]] = dependent$details$value[i]
      }


    }
    if(dependent$type == 'Grepl'){
      # In this case we use 'grepl()' to evalulate string patterns in the quantity.
      for(i in 1:nrow(dependent$details)){
        value[grepl(pattern = dependent$details$names[i],x = quantity)] = dependent$details$value[i]
      }
    }
    if(dependent$type == 'Number range (int)'){
      # In this case we use '%in%' and eval(parse) to check whether the numeric value is within a range of integers
      for(i in 1:nrow(dependent$details)){
        value[quantity %in% eval(parse(text = dependent$details$names[i]))] = dependent$details$value[i]
      }
    }
    if(dependent$type == 'Number range'){
      # In this case we use \leq and \geq.
      for(i in 1:nrow(dependent$details)){
        bounds = as.numeric(stringr::str_split(dependent$details$names[i], pattern = ':')[[1]])
        value[as.numeric(quantity) >= bounds[1] & as.numeric(quantity) <= bounds[2] ] = dependent$details$value[i]
      }
    }
    value = value * dependent$weight

    return(value)
  }))
  names(values_each_dependent) = names(dependents)

  ### return total value and breakdown of value.
  return(list(total_value = rowSums(values_each_dependent), value_breakdown = values_each_dependent))

  }
