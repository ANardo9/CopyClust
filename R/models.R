#Global Variables
utils::globalVariables(c("IntClust_Label", "IntClust_colors", "ID", "loc.end", "loc.start", "num.mark", "width",
                         "hg18_ranges", "hg19_ranges", "hg38_ranges", "Range", "Value"))
IntClust_colors = c("#E94D03","#7CB772","#B93377","#6EB8BB","#782D24","#F3E855","#364085","#E4AA2B","#E696E1","#6E3387")

#' CopyClust model integrative cluster prediction for breast cancer tumors
#'
#'@description
#' `CopyClust()` implements an XGBoost-based classifier trained on the copy number profiles of the METABRIC cohort to predict integrative cluster label based on copy number data alone. Integrative cluster prediction can be made using either a 10-class model approach
#' or 6-class with binary reclassification model approach. Scaling of features occurs prior to classification, therefore, accuracy may be impaired by data sets with small sample size.
#'
#' @param data_input A data frame with sample IDs as rows and the 478 un-scaled model features as columns. Can be the output from [CC_format()].
#' @param model_approach Parameter for model approach. If equal to `10C`, will implement 10-class model approach. If equal to `6C`, will implement 6-class model approach with binary reclassification. Default is `6C`, the 6-class with binary reclassification approach.
#' @returns A named numeric vector of predicted integrative cluster label the same length of number of samples provided with sample ID as row name.
#' @seealso [CC_format()] for a convenient way of formatting copy number data from `DNACopy` format into a structure usable by [CopyClust()].
#' @author Cameron C. Young
#' @examples
#' data("test_data")
#' results = CopyClust(test_data, model_approach = "6C")
#'
#' @export
#'
#' @importFrom stats predict
#' @importFrom tidyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom xgboost xgb.load.raw
#'

CopyClust = function(data_input, model_approach = "6C") {
  data_input = as.matrix(data_input)

  #Add error for incorrect data format
  if (dim(data_input)[2] != 478) {
    stop("Incorrect data_input format. Ensure data_input contains the 478 model features as columns and individual samples as rows.")
  }

  #10-Class Model
  if (model_approach == "10C") {
   data_input = scale(data_input)
   CopyClust_10_Class_Scale_model = xgb.load.raw(CopyClust_10_Class_Scale_model)
   prediction = as.data.frame(predict(CopyClust_10_Class_Scale_model, data_input) + 1)
   rownames(prediction) = rownames(data_input)
   colnames(prediction) = "IntClust_Label"

   return(prediction)
  }

  #6-Class Model with Binary Reclassification
  if (model_approach == "6C") {
    data_input = scale(data_input)
    CopyClust_6_Class_Scale_model = xgb.load.raw(CopyClust_6_Class_Scale_model)
    prediction = factor(predict(CopyClust_6_Class_Scale_model, data_input))
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
    data_15 = as.matrix(data_input[which(rownames(data_input) %in% rownames(samples_15)),])
    if(dim(data_15)[1] == 478 & dim(data_15)[2] == 1) {
      data_15 = t(data_15)
      rownames(data_15) = rownames(samples_15)
    }

    data_38 = as.matrix(data_input[which(rownames(data_input) %in% rownames(samples_38)),])
    if(dim(data_38)[1] == 478 & dim(data_38)[2] == 1) {
      data_38 = t(data_38)
      rownames(data_38) = rownames(samples_38)
    }

    data_47 = as.matrix(data_input[which(rownames(data_input) %in% rownames(samples_47)),])
    if(dim(data_47)[1] == 478 & dim(data_47)[2] == 1) {
      data_47 = t(data_47)
      rownames(data_47) = rownames(samples_47)
    }
    data_910 = as.matrix(data_input[which(rownames(data_input) %in% rownames(samples_910)),])
    if(dim(data_910)[1] == 478 & dim(data_910)[2] == 1) {
      data_910 = t(data_910)
      rownames(data_910) = rownames(samples_910)
    }

    #Implement Binary Models
    CopyClust_15_Binary_model = xgb.load.raw(CopyClust_15_Binary_model)
    samples_15_prediction = predict(CopyClust_15_Binary_model, data_15)
    samples_15_prediction = as.data.frame(ifelse(samples_15_prediction < 0.5, 1, 5))
    colnames(samples_15_prediction) = "IntClust_Label"
    rownames(samples_15_prediction) = rownames(data_15)

    CopyClust_38_Binary_model = xgb.load.raw(CopyClust_38_Binary_model)
    samples_38_prediction = predict(CopyClust_38_Binary_model, data_38)
    samples_38_prediction = as.data.frame(ifelse(samples_38_prediction < 0.5, 3, 8))
    colnames(samples_38_prediction) = "IntClust_Label"
    rownames(samples_38_prediction) = rownames(data_38)

    CopyClust_47_Binary_model = xgb.load.raw(CopyClust_47_Binary_model)
    samples_47_prediction = predict(CopyClust_47_Binary_model, data_47)
    samples_47_prediction = as.data.frame(ifelse(samples_47_prediction < 0.5, 4, 7))
    colnames(samples_47_prediction) = "IntClust_Label"
    rownames(samples_47_prediction) = rownames(data_47)

    CopyClust_910_Binary_model = xgb.load.raw(CopyClust_910_Binary_model)
    samples_910_prediction = predict(CopyClust_910_Binary_model, data_910)
    samples_910_prediction = as.data.frame(ifelse(samples_910_prediction < 0.5, 9, 10))
    colnames(samples_910_prediction) = "IntClust_Label"
    rownames(samples_910_prediction) = rownames(data_910)

    prediction = as.data.frame(rbind(samples_38_prediction, samples_47_prediction, samples_910_prediction, samples_15_prediction, other_samples))
    prediction = prediction[order(rownames(prediction)),,drop=FALSE]

    #Return prediction
    return(prediction)
  }
  else {
    stop("Incorrect input to 'model_approach' parameter. Must be equal to '6C' for 6-class with binary reclassification approach or '10C' for 10-class approach.")
  }
}

#' Format data for use by CopyClust function
#'
#' @description
#' Formats raw data from DNACopy format into 478 genomic range features required to run the [CopyClust()].
#' Reference genome (`hg18`, `h19`, or `hg38`) must be specified with the `reference_genome` parameter.
#'
#' @param data_input A data frame generated from the output of [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html). The data frame must contain six columns with the following column names: `ID`, `chrom`, `loc.start`, `loc.end`, `num.mark`, and `seg.mean`.
#' @param reference_genome Parameter for reference genome. Formats the genomic ranges used as features for [CopyClust()] to the appropriate reference genome. Valid inputs are `hg18`, `hg19`, and `hg38`. Default is `hg18`.
#' @param probes Parameter for number of probes to use to calculate feature values. Specified number of probes will be used to calculate values for model features selected from all available probes equally spaced across the genome. Default is `100000` probes. A greater number of probes decreases the processing speed. Providing a values greater than the number of available probes will results in a error.
#' @returns A data frame with sample IDs as rows and 478 model features as columns that can be used with the [CopyClust()] function.
#' @seealso [CopyClust()], [DNAcopy](https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)
#' @author Cameron C. Young
#' @examples
#' data("test_data_raw")
#' data_for_CopyClust = CC_format(test_data_raw, reference_genome = "hg19", probes = 100000)
#'
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
    message("hg18 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    message(paste("Number of samples to format: ", as.numeric(length(sample_ids)), sep = ""))

    #Expand data
    for(id_index in 1:length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Identify probe locations for use
      if(probes > sum(as.numeric(data_subset[,5]))) {
        stop(paste("Value of `probes` greater than number of available probes: ", sum(as.numeric(data_subset[,5])), sep = ""))
      }
      probe_locations = round(seq(from = 1, to = sum(as.numeric(data_subset[,5])), length.out = probes))

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = NA
      for (i in 1:dim(data_subset)[1]) {
        position_output_seq = seq(from = as.numeric(data_subset[i,3]), to = as.numeric(data_subset[i,4]), length.out = as.numeric(data_subset[i,5]))
        position_output = c(position_output, position_output_seq)
      }
      position_output = position_output[-1]

       #Add position to expanded_data
       expanded_data[,3] = position_output

       #Remove unused probes
       expanded_data = expanded_data[probe_locations,]

       #Which Range Loop
       range_output = numeric(length = dim(expanded_data)[1])
       for (i in 1:dim(expanded_data)[1]) {
         range_output[i] = hg18_ranges$range[max(which(hg18_ranges$start < as.numeric(expanded_data[i,3]) &
                                                         hg18_ranges$chrom == as.numeric(expanded_data[i,2])))]
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

        for (i in 1:dim(expanded_data)[1]) {
          feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
        }
         message(paste("Samples formatted: ", id_index, sep = ""))
     }
    return(feature_values)
  }
  else{

  #hg19
  if(reference_genome == "hg19") {
    message("hg19 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    message(paste("Number of samples to format: ", as.numeric(length(sample_ids)), sep = ""))

    #Expand data
    for(id_index in 1:length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Identify probe locations for use
      if(probes > sum(as.numeric(data_subset[,5]))) {
        stop(paste("Value of `probes` greater than number of available probes: ", sum(as.numeric(data_subset[,5])), sep = ""))
      }
      probe_locations = round(seq(from = 1, to = sum(as.numeric(data_subset[,5])), length.out = probes))

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = NA
      for (i in 1:dim(data_subset)[1]) {
        position_output_seq = seq(from = as.numeric(data_subset[i,3]), to = as.numeric(data_subset[i,4]), length.out = as.numeric(data_subset[i,5]))
        position_output = c(position_output, position_output_seq)
      }
      position_output = position_output[-1]

      #Add position to expanded_data
      expanded_data[,3] = position_output

      #Remove unused probes
      expanded_data = expanded_data[probe_locations,]

      #Which Range Loop
      range_output = numeric(length = dim(expanded_data)[1])
      for (i in 1:dim(expanded_data)[1]) {
        range_output[i] = hg19_ranges$range[max(which(hg19_ranges$start < as.numeric(expanded_data[i,3]) &
                                                        hg19_ranges$chrom == as.numeric(expanded_data[i,2])))]
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

      for (i in 1:dim(expanded_data)[1]) {
        feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
      }
      message(paste("Samples formatted: ", id_index, sep = ""))
    }
    return(feature_values)
  }
  else{

  #hg38
  if(reference_genome == "hg38") {
    message("hg38 Reference Genome")

    #Isolate samples IDs
    sample_ids = as.character(levels(factor(data_input$ID)))

    #Create matrix for feature values for all input samples
    feature_values = matrix(nrow = length(sample_ids), ncol = 478)
    rownames(feature_values) = sample_ids

    message(paste("Number of samples to format: ", as.numeric(length(sample_ids)), sep = ""))

    #Expand data
    for(id_index in 1:length(sample_ids)) {

      #Isolate single samples
      data_subset = data_input %>%
        mutate(width = loc.end - loc.start,
               separation = width / num.mark) %>%
        filter(ID == levels(factor(sample_ids))[id_index])
      data_subset = as.matrix(data_subset)

      #Identify probe locations for use
      if(probes > sum(as.numeric(data_subset[,5]))) {
        stop(paste("Value of `probes` greater than number of available probes: ", sum(as.numeric(data_subset[,5])), sep = ""))
      }
      probe_locations = round(seq(from = 1, to = sum(as.numeric(data_subset[,5])), length.out = probes))

      #Create expanded data matrix for isolated sample
      expanded_data = matrix(nrow = sum(as.numeric(data_subset[,5])), ncol = 5)
      colnames(expanded_data) = c("ID", "Chrom", "Position", "Value", "Range")

      expanded_data[,1] = as.character(rep(levels(factor(sample_ids))[id_index]), times = as.numeric(data_subset[,5]))
      expanded_data[,2] = as.numeric(rep(data_subset[,2], times = as.numeric(data_subset[,5])))
      expanded_data[,4] = as.numeric(rep(data_subset[,6], times = as.numeric(data_subset[,5])))

      #Position Loop
      position_output = NA
      for (i in 1:dim(data_subset)[1]) {
        position_output_seq = seq(from = as.numeric(data_subset[i,3]), to = as.numeric(data_subset[i,4]), length.out = as.numeric(data_subset[i,5]))
        position_output = c(position_output, position_output_seq)
      }
      position_output = position_output[-1]

      #Add position to expanded_data
      expanded_data[,3] = position_output

      #Remove unused probes
      expanded_data = expanded_data[probe_locations,]

      #Which Range Loop
      range_output = numeric(length = dim(expanded_data)[1])
      for (i in 1:dim(expanded_data)[1]) {
        range_output[i] = suppressWarnings(hg38_ranges$range[max(which(hg38_ranges$start < as.numeric(expanded_data[i,3]) &
                                                        hg38_ranges$chrom == as.numeric(expanded_data[i,2])))])
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

      for (i in 1:dim(expanded_data)[1]) {
        feature_values[as.numeric(id_index), as.numeric(expanded_data$Range[i])] = as.numeric(expanded_data$range_mean[i])
      }
      message(paste("Samples formatted: ", id_index, sep = ""))
    }
    return(feature_values)
  }
  else {
  stop("Incorrect entry to parameter 'reference_genome'.")
      }
    }
  }
}




