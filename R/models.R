#' CopyClust Model Prediction
#'
#' @param data_input A data frame with sample IDs as rows and 478 model features as columns.
#' @returns A numeric vector of predicted Integrative Cluster label.

CopyClust_Model_Prediction = function(data_input) {
   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_Function_v2, data_input)) + 1
   rownames(prediction) = rownames(data_input)
   colnames(prediction) = "IntClust Label"

   return(prediction)
}



