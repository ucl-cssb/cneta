#!/usr/bin/env Rscript

# suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

# This script includes functions to process the output of QDNASeq to prepare input to CNETML

## Transform the data in the list element to a data frame
transform.data <- function(mat, bins_data){
  df <- as.data.frame(t(mat))

  df <- cbind(df, bins_data)

  df <- df[, c("chromosome","start","end","copynumber","calls")]
  
  # assign bin ID by chr
  df = df %>% group_by(chromosome) %>% dplyr::mutate(bin = 1:n())
  return(df)
}


# p: patient ID
# cna_all: a list of CNA calls for all the samples, each list contains a data frame with at least five columns: "chromosome","start","end","copynumber","calls"
# bins: the boundary of bins of the CNA calls, with 3 columns: chromosome, start, end
# samps: all the samples of the patient, a data frame with 3 columns: "sample", "patient", "region"
# dir: output directory
write_cn <- function(p, cna_all, bins, samps, dir){
    nsamps <- nrow(samps)
    d <- vector("list", nsamps)
    # d[[i]] contains nbins x 5 ( bin, chr, start, end, cn )
    for(i in 1:nsamps){
        id <- samps[i,]$sample
        d[[i]] <- transform.data(cna_all[[id]], bins)
        #cat("variants:", unique(d[[i]]$calls), "\n")
        d[[i]]$id <- i
    }
    dd <- bind_rows(d)

    dd <- dd[, c("id", "chromosome", "bin", "calls")]
    dd$calls <- dd$calls + 2
    d.g <- data.frame(id=nsamps+1, chromosome=d[[1]]$chromosome, bin=d[[1]]$bin, calls=2)

    dd <- rbind(dd, d.g)

    ## write out data
    fout <- file.path(dir, paste("data-", p, "-cn.txt", sep=""))
    gz1 <- gzfile(paste(fout, ".gz", sep=""), "w")
    write.table(dd, gz1, quote=F, row.names=F, col.names=F, sep="\t")
    close(gz1)
}


# p: patient ID
# samps: all the samples of the patient, a data frame with 3 columns: sample, patient, region
# patient_info: a data frame with patient sample information, with columns including "PatientID and "sample", and additional columns when the sampling time information is available, "Age" and "SampleDate"
# incl_time: whether or not the sampling time information is included in patient_info
write_patient_info <- function(p, samps, patient_info, dir, incl_time = F){
  nsamps <- nrow(samps)

  ds <- data.frame(id=1:(nsamps+1), sample=c(as.character(samps$sample), "unaltered"))
  ds <- merge(ds, patient_info, by=c("sample"), all.x=T)
  # ds <- ds %>% select(id, sample, SampleDate, Age)
  ds[ds==""] = NA
  ds <- ds[order(ds$id),]
  
  sout <- file.path(dir, paste("data-", p, "-sample-ids.txt",sep=""))
  write.table(ds, sout, quote=F, row.names=F, col.names=T, sep="\t")

  if(incl_time){
    # sort the patient_info and calculate days
    ds %>% filter(is.na(SampleDate)) %>% nrow() -> n_nodate
    
    if(n_nodate == 1){
        ds$SampleDate <- as.Date(ds$SampleDate, format="%d/%m/%Y")
        mindate <- min(ds$SampleDate, na.rm=T)
        ds$Reldate <- difftime(ds$SampleDate, mindate, units="days") / 365
        # Get birth year at first, which is not in the excel table
        ds %>% filter(Age > 0) %>% head(1) -> tmp
        birth_year <- as.numeric(substring(tmp$SampleDate, 1, 4)) - as.numeric(tmp$Age)
        ds$Age <- as.numeric(substring(ds$SampleDate, 1, 4)) - birth_year  
    }else{
        ds$Reldate = 0
    }
  
    dt <- data.frame(x=c(1:nsamps), y=0)
    # sample, relative time, age
    dt <- ds %>% select(id, Reldate, Age)
    # remove the last row
    dt <- dt[1:(nrow(dt)-1),]

    tout <- file.path(dir, paste("data-", p, "-rel-times.txt",sep=""))
    write.table(dt, tout, quote=F, row.names=F, col.names=F, sep="\t")
  }
}


# cna_all: a list of CNA calls for all the samples, each list contains a data frame with at least five columns: "chromosome","start","end","copynumber","calls"
# bins: the boundary of bins of the CNA calls, with 3 columns: chromosome, start, end
# patient_info: a data frame with patient sample information, with columns including PatientID, sample, Age, SampleDate        
# patient_sample: a data frame with 3 columns: sample, patient, region
# patients: a vector of PatientIDs to process
# dir: output directory
write_patients <- function(cna_all, bins, patient_info, patient_sample, patients, dir, to_write_cn = T, to_write_info = T){ 
  for(p in patients){
    # p = patients[1]
    print(p)

    samps <- patient_sample[patient_sample$patient == p,]
    #cat("patient/sample:", p, as.character(samps$sample), samps$id,"\n")
    cat("patient/nsample:", p, length(as.character(samps$sample)), "\n")
    print(samps$sample)

    if(to_write_cn){
      write_cn(p, cna_all, bins, samps, dir)
    }
    
    if(to_write_info){
      write_patient_info(p, samps, patient_info, dir)
    }
  }
}



# TO ADD: sample command call
