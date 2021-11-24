layout: page
title: "Motif Analysis for Ecoacoustics"
permalink: https://ninascarpelli.github.io/MotifAnalysisComplete

# Motif Analysis for Ecoacoustics

The analysis developed here were motivated by the need of having a way to analyse multiple acoustic indices statistically. Currently, we have great visualisation tools but the step between those and the actual stats for ecological assessment was missing. The analysis showing here was accepted for publication and I will update the DOI and information on the paper as soon as I have them.

This code runs the motif analysis up to wavelet transformation with images and a way to check files for labelling the motifs. It also gives the samples to be labelled (30% of dataset). The second step is the file named 2_.... (Haven't done this yet whoops) and it performs the random forest classification. The second step should be run after at least 10% of motifs is labelled - ideally 30%.

If you're just after the script, you can find it [here](https://github.com/ninascarpelli/ninascarpelli.github.io/blob/a64347a74d83f3024efc019a2058b42346586737/1_MotifAnalysisComplete.R)

If you want the step-by-step explanation, here it is:

So, first here are the packages we will need for the analysis
```

library(tidyverse)
library(ggplot2)
library(wavelets)
library(dtwclust)
library(randomForest)
library(magick)
library(rJava)
```

Removing all the previous objects stored
```
rm(list = ls())

```


and some functions to make life easier when tasks are repetitive

1. First function is ```getDataPath``` that wasn't me who did it. This is so we can have the path to the parent folder that we will be working on without having to type the full path all the time. Best practice, according to Anthony (who is always trying to make me learn them) is to have absolute paths whenever possible. Also, create functions whenevee you have to copy and paste lines of code more than... 3 times?

```
getDataPath <- function (...) {
  return(file.path("C:",  ...))
}
```

2. ```list_myfiles``` is going to list all the files in one folder because most times we will be using for loops and doing the same task in different files.

```
list_myfiles <- function(step = NULL, search_pattern) {
  if (is.null(step)) {
    list.files(
      getDataPath(folder),
      pattern = search_pattern,
      recursive = T,
      full.names = T
    )
    
  } else {
    list.files(
      getDataPath(folder, step),
      pattern = search_pattern,
      recursive = T,
      full.names = T
    )
  }
  
}
```
3. ```create_mydir``` will create a directory inside the parent folder so we can store results
```
create_mydir <- function(current_step) {
  dir.create(getDataPath(folder, current_step))
}
```
4. ```ordering_files``` will be needed in step 2 when we build the time-series, i.e. putting files in order.

```
ordering_files <-
  function(original_dataframe,
           index_name) {
    if (is.null(index_name)) {
      c <- read.csv(original_dataframe) %>%
        with(., .[order(date_time, ResultMinute), ])
      new_dataframe <<- rbind(new_dataframe, c)
      write.csv(new_dataframe,
                getDataPath(folder,
                            step2,
                            paste("TS_", geo, "_", month, ".csv", sep = "")),
                row.names = F)
      
    } else {
      c <- read.csv(original_dataframe) %>%
        with(., .[order(date_time, ResultMinute), ]) %>%
        select(., all_of(index_name))
      new_dataframe <<- rbind(new_dataframe, c)
      write.table(
        new_dataframe,
        getDataPath(
          folder,
          step2,
          paste("TS_", geo,
                "_",
                month,
                "_",
                index_name,
                ".txt",
                sep = "")
        ),
        row.names = F,
        col.names = F
      )
      
      
    }
    
  }
```
4. ```HIME``` calls the HIME algorithm in powershell. This can be done separately by openning Powershell and running the command directly, I just put it here so it is easier for people without familiarity with Powershell.

```
HIME <- function(command) {
  system2('powershell', command)
}
```
5. This might be the most complicated function here, but nothing to worry about. ```iteration1``` will be used when cleaning up the motifs seeing there are lots of repetition. Again, the authors of the paper that did the algorithm and everything propose a different solution to cleaning them up, but this was the way I found to make it easier and straightforward when dealing with heaps of data as I was. This function is assigning the word "repeated" to motifs that overlap. 

```
iteration1 <- function(motif_results_df) {
  for (row in 1:nrow(motif_results_df)) {
    motif_results_df$overlap[row] = case_when(
      motif_results_df$Start[row + 1] <= motif_results_df$Start[row] &
        motif_results_df$End[row + 1] >= motif_results_df$End[row] ~ "repeated",
      motif_results_df$Start[row +
                               1] <= motif_results_df$End[row] &
        motif_results_df$End[row + 1] <= motif_results_df$End[row] ~ "repeated",
      motif_results_df$Start[row +
                               1] <= motif_results_df$End[row] &
        motif_results_df$End[row + 1] - motif_results_df$End[row] <= 30 ~ as.character(motif_results_df$Distance[row +
                                                                                                                   1]),
      motif_results_df$Start[row +
                               1] <= motif_results_df$End[row] &
        motif_results_df$Start[row + 1] - motif_results_df$Start[row] <= 30 ~ as.character(motif_results_df$Distance[row +
                                                                                                                       1]),
      motif_results_df$Start[row +
                               1] - motif_results_df$Start[row] <= 30 &
        motif_results_df$Start[row + 1] - motif_results_df$End[row] <= 30 ~ as.character(motif_results_df$Distance[row +
                                                                                                                     1]),
      TRUE ~ as.character(motif_results_df$Distance[row])
    )
  }
  motif_results_df <- motif_results_df
}
```

6. Now this function actually cleans up the motifs. ```remove_repeated``` runs the iteration 3 times so it removes all the possible repetitions and overlaps and remove motifs assigned as repeated.

```
remove_repeated <- function(motif_results_df) {
  result <- iteration1(motif_results_df) %>%
    filter(., .$overlap != "repeated") %>%
    with(., .[order(Start), ]) %>%
    iteration1(.) %>%
    filter(., .$overlap != "repeated") %>%
    with(., .[order(Start), ]) %>%
    iteration1(.) %>%
    filter(., .$overlap != "repeated") %>%
    with(., .[order(Start), ])
  
}
```

## Some more housekeeping before we start
Firstly name the parent folder where you have your indices and where we will be creating new folders to save our results
```
folder <- "AIndices"
```

### Steps and folders
It is good to have the intermediate steps saved so if something goes wrong you can fix it from where it stopped and doesn't have to run everything again.
The file structure is quite important bacause the way it is saved to avoid confusion. Therefore, files should be in a folder witn the site name, inside this folder should be the ones with the point names and then the files to be processed inside - following the pattern in the AP.exe output.

This analysis is made of 9 steps. Again, maybe that is not the best way of doing it, but this was the way I found best to have intermediate steps saved to be able to double check them. Also, the steps make sense in my head so it makes it straightforward to me. Maybe won't be the same for everyone so just change it however you think it is best for you.

The steps are:
```
step1 <- "1_IndicesToTs"
step2 <- "2_MonthlyAI"
step3 <- "3_HIME"
step4 <- "4_ProcessingRes"
step5 <- "5_CleaningUp"
step6 <- "6_CompleteMotif"
step7 <- "7_CropSpectrogram"
step8 <- "8_FeatureExtraction"
```

And then we create directories for all those steps
```
create_mydir(step1)
create_mydir(step2)
create_mydir(step3)
create_mydir(step4)
create_mydir(step5)
create_mydir(step6)
create_mydir(step7)
create_mydir(step8)
create_mydir("Figures")
```

## Nice, we are finally ready to start
### Step 1: ```1_IndicesToTs```
Here we will ordenate the indices minutes as a time-series

Listing the files in the directory:
```
files <- list_myfiles(search_pattern = ".Indices.csv")
```

Creating dirs and building TS 
First create and unique identifier to each minute using the file name - which ideally has date and time embedded - and then scale the indices values - if necessary adjust the columns in line 15 - scale only columns with indices

This create one file with all the 3 indices (we will use the FCS ones) + all variables needed for subsequence analysis. It also saves the sites, points and date_time attributes so we use it in the next step

```
site_id <- NULL
point_id <- NULL
date_time_id <- NULL
month_id <- NULL
geo_id <- NULL

for (file in files) {
  t <- strsplit(file, split = "/")
  u <- strsplit(t[[1]][6], split = "_")
  d <- substr(u[[1]][1], start = 1, stop = 6)
  
  geo_id <- c(geo_id, paste(t[[1]][3], t[[1]][4], sep = "_"))
  site_id <- c(site_id, t[[1]][3])
  point_id <- c(point_id, t[[1]][4])
  date_time_id <- c(date_time_id, u[[1]][1])
  month_id <- c(month_id, d[1])
  
  site_id <- unique(site_id)
  point_id <- unique(point_id)
  date_time_id <- unique(date_time_id)
  month_id <- unique(month_id)
  geo_id <- unique(geo_id)
  
  
  if (dir.exists(t[[1]][3])) {
    cat("The folder already exists")
    
  } else {
    dir.create(getDataPath(folder, step1, t[[1]][3]))
    
  }
  
  if (dir.exists(t[[1]][4])) {
    cat("The folder already exists")
    
  } else {
    dir.create(getDataPath(folder, step1, t[[1]][3], t[[1]][4]))
    
  }
  
  
  read.csv(file) %>%
    mutate(., FID = paste(basename(file), ResultMinute, sep = "_")) %>%
    separate(
      .,
      FID,
      into = c("date_time", "useless", "useless1", "useless2"),
      sep = "_",
      remove = F
    ) %>%
    separate(
      .,
      date_time,
      into = c("date", "time"),
      sep = "T",
      remove = F
    ) %>%
    mutate(., site = t[[1]][3]) %>%
    mutate(., point = t[[1]][4]) %>%
    mutate(., filepath = file) %>%
    mutate(., month = d[1]) %>%
    select(
      .,
      AcousticComplexity,
      EventsPerSecond,
      TemporalEntropy,
      FileName,
      date_time,
      month,
      date,
      time,
      ResultMinute,
      FID,
      site,
      point,
      filepath
    ) %>%
    mutate_at(vars(1:3), scale) %>%
    with(., .[order(as.numeric(date), as.numeric(time), ResultMinute),]) %>%
    write.csv(., getDataPath(folder,
                             step1,
                             t[[1]][3],
                             t[[1]][4],
                             paste(t[[1]][3], "_", t[[1]][4], "_", u[[1]][1], ".csv", sep = "")), row.names = F)
}

```
### Step 2: ```2_MonthlyAI```

This analysis works for resolution of 24 hours up to one month. We did it this way because the motifs try to find repetitions in the time-series and we expect some ecological patterns to be repeated across different time-scales. However we wanted to keep monthly data separated respecting the natural cycles that can occur with temperature, rainfall and many other changes. This was done so it would be -theoretically- easier to explain and understand patterns in a monthly resolution than across a year or years of recordings.

The output here should be one df per month/index = 3 dfs per month.

```
new_dataframe <- NULL

for (geo in geo_id) {
  for (month in month_id) {
    files <-
      list_myfiles(step1,
                   search_pattern = paste(geo, "_", month, sep = ""))
    
    for (file in files) {
      ordering_files(file, index_name = NULL)
    }
    new_dataframe <- NULL
    
    for (file in files) {
      ordering_files(file, index_name = "AcousticComplexity")
    }
    
    new_dataframe <- NULL
    for (file in files) {
      ordering_files(file, index_name = "TemporalEntropy")
    }
    
    new_dataframe <- NULL
    for (file in files) {
      ordering_files(file, index_name = "EventsPerSecond")
    }
    
    new_dataframe <- NULL
    
  }
}
```

### Step 3: ```3_HIME```

Set the directory containing the files
```
input_directory <- getDataPath(folder, step2)
```
The directory to store the results
```
output_directory <- getDataPath(folder, step3)
```
Build the command

```
HIME_command <-
  paste("'", input_directory, "/", "HIME_release.jar", "'", sep  = "")
```

Command to set directory in PowerShell

```
command1 <- paste("cd ", input_directory, sep = "")
```
Listing files
```
files <- list_myfiles(step2, search_pattern = glob2rx("TS_*.txt"))
```
Iterate through each file, prepare the command and execute in the for loop
```
for (file in files) {
  pws_input_file <-
    paste("'", input_directory, "/", basename(file), "'", sep = "")
  pws_output_file <-
    paste("'",
          output_directory,
          "/",
          gsub("TS_", "res", basename(file)),
          "'",
          sep = "")
  
  # prepare command
  
  command2 <-
    paste("java -jar", HIME_command, pws_input_file, "4 32 > tmp.log")
  command3 <-
    paste("select-string Motif tmp.log >", pws_output_file)
  
  # finally, execute the command
  HIME(command1)
  HIME(command2)
  HIME(command3)
}

```

### Step4: ```4_ProcessingRes```

After motif analysis, this will process .txt file from motif - changin encoding and getting rid of columns and info we don't need

```
files <- list_myfiles(step3, search_pattern = glob2rx("res*.txt"))

for (file in files) {
  read.table(
    file,
    sep = " ",
    blank.lines.skip = T,
    fileEncoding = "UTF-16"
  ) %>%
    select(., 2:7) %>%
    write.table(getDataPath(folder, step4, gsub("res", "", basename(file))),
                row.names = F,
                col.names = F)
}
```

### Step5: ```5_CleaningUp```

This step will create the plots with the time series and motifs + clean them up removing the repetitions I mentioned before

Firstly, naming the indices for the plots
```
indices <-
  c("AcousticComplexity", "EventsPerSecond", "TemporalEntropy")
```
Plotting TS for all indices
```
for (geo in geo_id) {
  for (month in month_id) {
    
    skip_to_next <- FALSE
    
    tryCatch({
    complete_ts <-
      read.csv(getDataPath(folder,
                           step2,
                           paste("TS_", geo, "_", month, ".csv", sep = ""))) %>%
      mutate(., position = seq_len(nrow(.)))
    
    },
    
    error = function(e) {
      skip_to_next <<- TRUE
    })
    
    if (skip_to_next) {
      next
      
    }
    
    plot_ts <-
      select(complete_ts, position, all_of(indices)) %>%
      pivot_longer(.,
                   cols = 2:4,
                   names_to = "index",
                   values_to = "value")
    
    ggplot(plot_ts, aes(x = position, y = value)) +
      geom_line() +
      facet_wrap(. ~ index) +
      theme_classic() +
      theme(axis.text.x = element_blank()) +
      ggsave(getDataPath(
        folder,
        "Figures",
        paste(geo, month, "indicespertime.jpg", sep = "_")
      ))
        
    for (index in indices) {
      complete_inter <-
        select(complete_ts, all_of(index), 4:ncol(complete_ts)) %>%
        mutate(., motif = NA) %>%
        mutate(., distance = NA) %>%
        mutate(., length = NA) %>%
        mutate(., reference = "0_ts") %>%
        mutate(., id = 0) %>%
        rename(., Index = index)
      
      
      
      # Processing results
      
      motif_results <-
        read.table(getDataPath(
          folder,
          step4,
          paste(geo, "_", month, "_", index, ".txt", sep = "")
        )) %>%
        rename(
          .,
          FirstInstance_Start = V1,
          FirstInstance_End = V2,
          SecondInstance_Start = V3,
          SecondInstance_End = V4,
          Length = V5,
          Distance = V6
        ) %>%
        mutate(., id = 1:as.numeric(count(.))) %>%
        filter(., Distance <= 5) %>%
        select(., id, everything()) %>%
        pivot_longer(.,
                     cols = 2:5,
                     names_to = "Instance",
                     values_to = "position") %>%
        mutate(.,
               Instance = gsub(
                 pattern = "FirstInstance",
                 replacement = "motif",
                 x = Instance
               )) %>%
        mutate(.,
               Instance = gsub(
                 pattern = "SecondInstance",
                 replacement = "match",
                 x = Instance
               )) %>%
        separate(.,
                 Instance,
                 into = c("instance", "moment"),
                 sep = "_") %>%
        pivot_wider(., names_from = moment, values_from = position) %>%
        mutate(., instance = paste(id, instance, sep = "_")) %>%
        with(., .[order(Start), ]) %>%
        mutate(., overlap = NA) %>%
        remove_repeated(.)
        
      
      
      for (row in 1:nrow(complete_inter)) {
        

        skip_to_next <- FALSE
        
        tryCatch({
          complete_inter[motif_results$Start[row]:motif_results$End[row], c("motif", "distance", "length")] <-
            motif_results[row, c("instance", "Distance", "Length")] 
        },
        
        error = function(e) {
          skip_to_next <<- TRUE
        })
        
        if (skip_to_next) {
          next
        
        }
      }
      
      complete_inter <-  group_by(complete_inter, motif) %>% 
        add_count(.) %>% 
        filter(n>=30)
      
      
      write.csv(complete_inter,
                getDataPath(
                  folder,
                  step5,
                  paste(geo, month, index, "motif.csv", sep = "_")
                ),
                row.names = F)
      
      plot_ts <-
        select(complete_inter, reference, position, Index, date, time) %>%
        separate(.,
                 reference,
                 into = c("number", "what"),
                 remove = F)
      
      plot_motif <-
        select(complete_inter, motif, position, Index) %>%
        rename(., reference = motif) %>%
        filter(reference != "NA") %>%
        separate(.,
                 reference,
                 into = c("number", "what"),
                 remove = F)
      
      
      #6 in BNE = +10
      line_intercept1 <- filter(plot_ts, grepl("160000*", time)) %>%
        .[!duplicated(.$date), ] %>%
        mutate(time = 060000) %>%
        select(time, position)
      
      #18 in BNE = +10
      line_intercept2 <- filter(plot_ts, grepl("040000*", time)) %>%
        .[!duplicated(.$date), ] %>%
        mutate(time = 180000) %>%
        select(time, position)
      
      ggplot(plot_ts, aes(x = position, y = Index)) +
        geom_line(aes(colour = what, linetype = what), colour = "grey") +
        geom_vline(xintercept = line_intercept1$position, linetype = "dotted") +
        geom_text(data = line_intercept1,
                  aes(label = time, y = 10, size = 1),
                  check_overlap = T) +
        geom_vline(xintercept = line_intercept2$position, linetype = "dotted") +
        geom_text(data = line_intercept2,
                  aes(label = time, y = 10, size = 1),
                  check_overlap = T) +
        scale_linetype_manual(values = "dotted") +
        geom_line(data = plot_motif, aes(
          x = position,
          y = Index,
          colour = reference
        )) +
        scale_color_manual(values = c(replicate(nrow(
          motif_results
        ), "#2ca25f"))) +
        theme_classic() +
        labs(title = paste(index, sep = " ")) +
        theme(
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"
        ) +
        ggsave(getDataPath(
          folder,
          "Figures",
          paste(geo, "_", month, "_", index, "ts_motifs.jpg", sep = "")
        ))
      
    }
    
  }
}

```
### Step6: ```6_CompleteMotif```
### Step7: ```7_CropSpectrogram```
### Step8: ```8_FeatureExtraction```

Steps 6, 7 and 8 go together and they create a file with all the motifs; crop the spectrograms so we can use them in the labelling process and then do the wavelet feature extraction that will be used in the random forest algorithm

```
for (geo in geo_id) {
  for (month in month_id) {
    
    motif_complete <- NULL
    file_result <- NULL
    img_prep <- NULL
    ts_data <- NULL
    ts_list <- NULL
    
    
    motif_complete <- data.frame(
      position =	integer(),
      index_value = numeric(),
      FileName = factor(),
      date = integer(),
      time = integer(),
      ResultMinute = integer(),
      FID = character(),
      distance = numeric(),
      length = integer(),
      reference = character(),
      id	= character(),
      fid_what = factor(),
      site = factor(),
      point = factor(),
      date_time = factor()
    )
    
    files <-
      list_myfiles(step5,
                   search_pattern = glob2rx(paste(geo, "*motif.csv", sep = "_")))
    
    for (file in files) {
      file_result <- read.csv(file) %>%
        dplyr::filter(., motif != is.na(T)) %>%
        dplyr::rename(., fid_what = motif) %>%
        dplyr::rename(., index_value = Index) %>%
        dplyr::mutate(.,
                      id = paste(
                        basename(file) %>%
                          gsub(pattern = "*_motif.csv", replacement = "") %>%
                          gsub(pattern = "TemporalEntropy", replacement = "ENT") %>%
                          gsub(pattern = "AcousticComplexity", replacement = "ACI") %>%
                          gsub(pattern = "EventsPerSecond", replacement = "EVN"),
                        fid_what,
                        sep = "_"
                      )) %>%
        dplyr::select(
          .,
          position,
          index_value,
          FileName,
          date,
          time,
          ResultMinute,
          distance,
          length,
          reference,
          id,
          fid_what,
          site,
          point,
          date_time
        )
      motif_complete <- rbind(motif_complete, file_result) %>% 
        filter(.$index_value != "NA")
      
img_prep <- separate(motif_complete,
        id,
        into = c("site", "point", "month", "index_name", "motif_number", "what"),
        remove = F
      ) %>%
        group_by(., id) %>%
        mutate(., new_position = order(order(position))) %>%
        ungroup(.) %>%
        select(everything(), -c(position)) %>%
        group_by(id) %>%
        filter(ResultMinute == min(ResultMinute))

dir.create(getDataPath(folder, step7, unique(img_prep$site)))
dir.create(getDataPath(folder, step7,unique(img_prep$site), unique(img_prep$point)))
      
      for (row in 1:nrow(img_prep)) {
        
        #list.files(getDataPath(folder, img_prep$site[row], img_prep$point[row]), pattern = glob2rx(paste(img_prep$month[row], "*", img_prep$index_name[row], ".png", sep = "")), recursive = T, full.names = T) %>%
       image_read(getDataPath(folder, img_prep$site[row], img_prep$point[row], paste(img_prep$FileName[row], "__", img_prep$index_name[row], ".png", sep = ""))) %>%
          image_crop(
            .,
            geometry_area(
              height = 256,
              width = img_prep$length[row] - (1 - img_prep$ResultMinute[row]),
              y_off = 20,
              x_off = img_prep$ResultMinute[row]
            )
          ) %>%
          image_write(., getDataPath(
            folder,
            step7,
            img_prep$site[row],
            img_prep$point[row],
            paste(img_prep$id[row],
              ".png",
              sep = ""
            )
          ))
      }
    }
    
    ts_data <- select(motif_complete, index_value, position, id) %>%
      group_by(., id) %>%
      mutate(., new_position = order(order(position))) %>%
      ungroup(.) %>%
      select(., everything(), -position) %>%
      pivot_wider(., names_from = new_position, values_from = index_value) %>%
      as.data.frame(.)

    
    
    rownames(ts_data) <- ts_data$id
    ts_data <- ts_data[,2:length(ts_data)]
    
    
    ts_list <- tslist(ts_data) %>%
      map(., na.omit)
    
    
    wtData <- NULL
    
    
    for (i in ts_list) {
      
      wt <- dwt(i, filter="haar", boundary= "periodic")
      
      un <- as.data.frame(t(unlist(c(wt@W,wt@V[[wt@level]]))))
      
      wtData <- plyr::rbind.fill(wtData, un)
      
    }
    
    wtData <- na.roughfix(wtData)
    
    wtData$id <- rownames(ts_data)
    
    wtData <- mutate(wtData, class = NA) %>% 
      mutate(., component = NA) %>% 
      select(., id, class, component, everything())
    
    
    samples <- sample(wtData$id, size = ceiling(nrow(wtData)*0.30), replace = F)
    
    write.csv(wtData, getDataPath(folder, step8, paste(geo, "_", month, "_wavelet.csv", sep = "")), row.names = F)
    write.csv(samples, getDataPath(folder, step8, paste(geo, "_", month, "_LabelsSample.csv", sep = "")), row.names = F)
    
    write.csv(motif_complete,
              getDataPath(
                folder,
                step6,
                paste(geo, month, "motif_complete.csv", sep = "_")
              ),
              row.names = F)

    
  }
}
```

Finally now you can label the motifs to train and run the random forest algorithm here.

