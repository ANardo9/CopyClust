#' CopyClust Model Prediction
#'
#' @param data_input A data frame with sample IDs as rows and 478 model features as columns.
#' @param model_approach If equal to "10C", implement 10-class model approach. Otherwise, will implement 6-class model approach with binary reclassification.
#' @returns A numeric vector of predicted Integrative Cluster label according to model approach.

CopyClust = function(data_input, model_approach = "10C") {
  if (model_approach == "10C") {

   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_Function_v2, data_input)) + 1
   rownames(prediction) = rownames(data_input)
   colnames(prediction) = "IntClust_Label"

   return(prediction)
  }
  else {
    prediction = factor(predict(CopyClust_6_Class_Scale_Function_v2, data_input))
    levels(prediction) = c("1/5", "2", "3/8", "4/7", "6", "9/10")
    prediction = as.data.frame(prediction)
    colnames(prediction) = "IntClust_Label"
    rownames(prediction) = rownames(data_input)

    #Filter Binary Groups
    samples_15 = prediction %>% filter(IntClust_Label == "1/5")
    samples_38 = prediction %>% filter(IntClust_Label == "3/8")
    samples_47 = prediction %>% filter(IntClust_Label == "4/7")
    samples_910 = prediction %>% filter(IntClust_Label == "9/10")
    other_samples = prediction %>% filter(IntClust_Label == "2" | IntClust_Label == "6")

    #Filter Binary Data
    data_15 = data_input[which(rownames(data_input) %in% rownames(samples_15)),]
    data_38 = data_input[which(rownames(data_input) %in% rownames(samples_38)),]
    data_47 = data_input[which(rownames(data_input) %in% rownames(samples_47)),]
    data_910 = data_input[which(rownames(data_input) %in% rownames(samples_910)),]

    #Implement Binary Models
    samples_15_prediction = predict(CopyClust_15_Binary_Scale_Function_v2, data_15)
    samples_15_prediction = as.data.frame(ifelse(samples_15_prediction < 0.5, 1, 5))
    colnames(samples_15_prediction) = "IntClust_Label"
    rownames(samples_15_prediction) = rownames(data_15)

    samples_38_prediction = predict(CopyClust_38_Binary_Scale_Function_v2, data_38)
    samples_38_prediction = as.data.frame(ifelse(samples_38_prediction < 0.5, 3, 8))
    colnames(samples_38_prediction) = "IntClust_Label"
    rownames(samples_38_prediction) = rownames(data_38)

    samples_47_prediction = predict(CopyClust_47_Binary_Scale_Function_v2, data_47)
    samples_47_prediction = as.data.frame(ifelse(samples_47_prediction < 0.5, 4, 7))
    colnames(samples_47_prediction) = "IntClust_Label"
    rownames(samples_47_prediction) = rownames(data_47)

    samples_910_prediction = predict(CopyClust_910_Binary_Scale_Function_v2, data_910)
    samples_910_prediction = as.data.frame(ifelse(samples_910_prediction < 0.5, 9, 10))
    colnames(samples_910_prediction) = "IntClust_Label"
    rownames(samples_910_prediction) = rownames(data_910)

    prediction = as.data.frame(rbind(samples_38_prediction, samples_47_prediction, samples_910_prediction, samples_15_prediction, other_samples))
    prediction = prediction[order(rownames(prediction)),,drop=FALSE]

    #Return prediction
    return(prediction)
    }
}









