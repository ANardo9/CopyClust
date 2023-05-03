#' CopyClust Model Prediction
#'
#' @param data_input A data frame with sample IDs as rows and 478 model features as columns.
#' @param model_approach If TRUE, implement 10-class model approach. If FALSE, implement 6-class model approach with binary reclassification.
#' @returns A numeric vector of predicted Integrative Cluster label according to model approach.

CopyClust_Model_Prediction = function(data_input, model_approach = TRUE) {
  if (model_approach == TRUE) {

   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_Function_v2, data_input)) + 1
   rownames(prediction) = rownames(data_input)
   colnames(prediction) = "IntClust Label"

   return(prediction)
  }
  else {
    print("6-Class Model Approach")
  }
}



