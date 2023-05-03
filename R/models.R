



CopyClust_10_Class_Scale_Model_Prediction = function(data_input) {
   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_Function_v2, data_input))
   return(prediction)
}



