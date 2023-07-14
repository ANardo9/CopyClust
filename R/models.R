#Global Variables
utils::globalVariables(c("IntClust_Label", "IntClust_colors", "ID", "loc.end", "loc.start", "num.mark", "width",
                         "hg18_ranges", "hg19_ranges", "hg38_ranges", "Range", "Value"))
IntClust_colors = c("#E94D03","#7CB772","#B93377","#6EB8BB","#782D24","#F3E855","#364085","#E4AA2B","#E696E1","#6E3387")

#' CopyClust Model Prediction
#'
#'@description
#'Implements XGBoost models developed using METABRIC cohort for IntClust classification based on copy number data alone.
#'
#' @param data_input A data frame with sample IDs as rows and 478 model features as columns.
#' @param model_approach Parameter for model approach, default is "6C": 6-Class with Binary Reclassification. If equal to "10C", implement 10-class model approach. If equal to "6C", will implement 6-class model approach with binary reclassification.
#' @returns A numeric vector of predicted Integrative Cluster label according to model approach.
#' @export
#'
#' @importFrom stats predict
#' @importFrom tidyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr group_by

CopyClust = function(data_input, model_approach = "6C") {
  #Add error for incorrect data format
  if (dim(data_input)[2] != 478) {
    stop("Incorrect data_input format. Ensure data_input contains 478 model features as columns and individual samples as rows.")
  }

  #10-Class Model
  if (model_approach == "10C") {
   data_input = scale(data_input)
   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_Function_v2, data_input)) + 1
   rownames(prediction) = rownames(data_input)
   colnames(prediction) = "IntClust_Label"

   return(prediction)
  }

  #6-Class Model with Binary Reclassification
  if (model_approach == "6C") {
    data_input = scale(data_input)
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
  else {
    stop("Incorrect input to 'model_approach' parameter.")
  }
}

#' Format Data for CopyClust Function
#'
#' @description
#' Formats raw data in DNACopy format into 478 genomic range features required to run the CopyClust() function.
#' Reference genome (hg18, h19, or hg38) must be specified.
#'
#' @param data_input A data frame representing the output of DNAcopy. Six columns: "ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"
#' @param reference_genome Formats the genomic ranges to the appropriate reference genome. Valid inputs are "hg18", "hg19", and "hg38". Default is "hg18".
#' @param probes Number of probes to utilize. Default is 100,000. A greater number of probes decreases the speed.
#' @returns A data frame with sample IDs as rows and 478 model features as columns that can be used with the CopyClust function.
#' @export
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup


CC_format = function(data_input, reference_genome = "hg18", probes = 100000) {
  #incorrect probes input
  if(!is.numeric(probes)){
    stop("Non-numeric entry to parameter 'probes'.")
  }

  #hg18
  if(reference_genome == "hg18") {
    print("hg18 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    #Expand data
    id_index = 1
    while(id_index <= length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = numeric(length = sum(as.numeric(data_subset[,5])))
      i = 1
      k = 1
      while (i <= dim(data_subset)[1]) {
        j = 1
          while (j <= as.numeric(data_subset[i,5])) {
            position_output[k] = as.numeric(data_subset[i,3]) + as.numeric(data_subset[i,8]) * (j - 1)
            j = j + 1
            k = k + 1
          }
        i = i + 1
      }

       #Add position to expanded_data
       expanded_data[,3] = as.numeric(round(position_output))

       #Which Range Loop
       range_output = numeric(length = sum(as.numeric(data_subset[,5])))
       i = 1
       while (i <= dim(expanded_data)[1]) {
         range_output[i] = hg18_ranges$range[max(which(hg18_ranges$start < as.numeric(expanded_data[i,3]) &
                                                         hg18_ranges$chrom == as.numeric(expanded_data[i,2])))]
         i = i + 1
       }

       #Add range to expanded_data
       expanded_data[,5] = as.numeric(round(range_output))

       #Calculate mean value by range
       expanded_data = as.data.frame(expanded_data)
       expanded_data$Range = as.numeric(expanded_data$Range)
       expanded_data$Value = as.numeric(expanded_data$Value)

       expanded_data = expanded_data %>%
         group_by(Range) %>%
         mutate(range_mean = mean(Value)) %>%
         ungroup()
       expanded_data = expanded_data[!duplicated(expanded_data$Range),]

        i = 1
        while (i <= dim(expanded_data)[1]) {
          feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
          i = i + 1
        }

        print(paste("Samples formatted: ", id_index, sep = ""))
        id_index = id_index + 1
    }
    return(feature_values)
  }
  else{

  #hg19
  if(reference_genome == "hg19") {
    print("hg19 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    #Expand data
    id_index = 1
    while(id_index <= length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = numeric(length = sum(as.numeric(data_subset[,5])))
      i = 1
      k = 1
      while (i <= dim(data_subset)[1]) {
        j = 1
        while (j <= as.numeric(data_subset[i,5])) {
          position_output[k] = as.numeric(data_subset[i,3]) + as.numeric(data_subset[i,8]) * (j - 1)
          j = j + 1
          k = k + 1
        }
        i = i + 1
      }

      #Add position to expanded_data
      expanded_data[,3] = as.numeric(round(position_output))

      #Which Range Loop
      range_output = numeric(length = sum(as.numeric(data_subset[,5])))
      i = 1
      while (i <= dim(expanded_data)[1]) {
        range_output[i] = hg19_ranges$range[max(which(hg19_ranges$start < as.numeric(expanded_data[i,3]) &
                                                        hg19_ranges$chrom == as.numeric(expanded_data[i,2])))]
        i = i + 1
      }

      #Add range to expanded_data
      expanded_data[,5] = as.numeric(round(range_output))

      #Calculate mean value by range
      expanded_data = as.data.frame(expanded_data)
      expanded_data$Range = as.numeric(expanded_data$Range)
      expanded_data$Value = as.numeric(expanded_data$Value)

      expanded_data = expanded_data %>%
        group_by(Range) %>%
        mutate(range_mean = mean(Value)) %>%
        ungroup()
      expanded_data = expanded_data[!duplicated(expanded_data$Range),]

      i = 1
      while (i <= dim(expanded_data)[1]) {
        feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
        i = i + 1
      }

      print(paste("Samples formatted: ", id_index, sep = ""))
      id_index = id_index + 1
    }
    return(feature_values)
  }
  else{

  #hg38
  if(reference_genome == "hg38") {
    print("hg38 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    #Expand data
    id_index = 1
    while(id_index <= length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = numeric(length = sum(as.numeric(data_subset[,5])))
      i = 1
      k = 1
      while (i <= dim(data_subset)[1]) {
        j = 1
        while (j <= as.numeric(data_subset[i,5])) {
          position_output[k] = as.numeric(data_subset[i,3]) + as.numeric(data_subset[i,8]) * (j - 1)
          j = j + 1
          k = k + 1
        }
        i = i + 1
      }

      #Add position to expanded_data
      expanded_data[,3] = as.numeric(round(position_output))

      #Which Range Loop
      range_output = numeric(length = sum(as.numeric(data_subset[,5])))
      i = 1
      while (i <= dim(expanded_data)[1]) {
        range_output[i] = hg38_ranges$range[max(which(hg38_ranges$start < as.numeric(expanded_data[i,3]) &
                                                        hg38_ranges$chrom == as.numeric(expanded_data[i,2])))]
        i = i + 1
      }

      #Add range to expanded_data
      expanded_data[,5] = as.numeric(round(range_output))

      #Calculate mean value by range
      expanded_data = as.data.frame(expanded_data)
      expanded_data$Range = as.numeric(expanded_data$Range)
      expanded_data$Value = as.numeric(expanded_data$Value)

      expanded_data = expanded_data %>%
        group_by(Range) %>%
        mutate(range_mean = mean(Value)) %>%
        ungroup()
      expanded_data = expanded_data[!duplicated(expanded_data$Range),]

      i = 1
      while (i <= dim(expanded_data)[1]) {
        feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
        i = i + 1
      }

      print(paste("Samples formatted: ", id_index, sep = ""))
      id_index = id_index + 1
    }
    return(feature_values)
  }
  else {
  stop("Incorrect entry to parameter 'reference_genome'.")
      }
    }
  }
}




