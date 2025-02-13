library(reticulate)
library(gtools)
library(Raff)
library(hadron)
library(stringr)

# Define Gamma_0 and Identity
g0 <- matrix(data = c(0,0,-1,0,0,0,0,-1,-1,0,0,0,0,-1,0,0), nrow = 4, ncol = 4, byrow = TRUE)
id <- matrix(data = c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), nrow = 4, ncol = 4, byrow = TRUE)

# Get the working directory
pwd <- dirname(sys.frame(1)$ofile)

# Source python files for NjjN calculations
source_python(paste(pwd,"/python/npv_classes.py",sep=""))
source_python(paste(pwd,"/python/npv_calc.py",sep=""))

# Write momentum in the key format
format_momentum <- function(momentum) {
  # Get the indices of negative momenta
  negative_momentum_ii <- which(momentum<0)
  # Absolute value of all 3-momentum
  abs_momentum <- abs(momentum)
  # Write in desired format i.e. "00", "01"
  for (i in 1:3) {
    momentum[i] <- sprintf("%02d",abs_momentum[i])
  }
  # Insert negative sign (-) where needed
  for (i in 1:length(negative_momentum_ii)) {
    momentum[negative_momentum_ii[i]] <- paste("-",momentum[negative_momentum_ii[i]],sep="")
  }
  return(momentum)
}

# Generate Key for Nucleon-Nucleon Correlator
generate_key_NN <- function(nucleon_type, source_coord, Gi = "Cg5", Gf = "Cg5", nucleon_momentum) {
  # Set the propagator ordering based on nucleon type
  if (nucleon_type == 'p'){ 
    k1 <- 'udu'
  }  else if (nucleon_type == 'n') { 
    k1 <- 'dud'
  }  else {
    print("nucleon_type must be 'n' or 'p'")
    return(NA)
  }
  # Write the keys for T1 and T2
  key1 <- paste("/N-N/",k1,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/Gi_",Gi,"/Gf_",Gf,"/t1/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
  key2 <- paste("/N-N/",k1,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/Gi_",Gi,"/Gf_",Gf,"/t2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
  key <- c(key1,key2)
  return(key)
}

# Get source coordinate from the file name
get_source_coord <- function(file_name) {
  # Split the filename into characters
  a <- strsplit(file_name,"")[[1]]
  # Find out the position of the second t
  # t <- which(a=="t")[2]+1
  t <- which(a=="t")[3]+1
  #set the t coordinate by reading the next position
  source_coord_t <- a[t]
  # If the position following is not "x" then the t coordinate is 2 digits, read and concatinate
  if (a[t+1]!="x") source_coord_t <- paste(source_coord_t,a[t+1],sep="")
  # similarly for "x", "y" and "z"
  x <- which(a=="x")[1]+1
  source_coord_x <- a[x]
  if (a[x+1]!="y") source_coord_x <- paste(source_coord_x,a[x+1],sep="")
  y <- which(a=="y")[1]+1
  source_coord_y <- a[y]
  if (a[y+1]!="z") source_coord_y <- paste(source_coord_y,a[y+1],sep="")
  z <- which(a=="z")[1]+1
  source_coord_z <- a[z]
  if (a[z+1]!=".") source_coord_z <- paste(source_coord_z,a[z+1],sep="")
  source_coord <- as.integer(c(source_coord_t,source_coord_x,source_coord_y,source_coord_z))
  return(source_coord)
}

# Calculate Nucleon-Nucleon correlator
NN_corr <- function(file_name, file_dir = "/home/aniket/projects/npv_data/cA211a.30.32/", nucleon_type, Gi = "Cg5", Gf = "Cg5", 
                    nucleon_momentum = c(0,0,0), LT = 64, parity = +1, read_from_disk = FALSE, write_to_disk = FALSE, read_write_dir = NA) {
  # Split the file name to get the extension
  f <- strsplit(file_name, split = "")[[1]]
  f_ext <- paste(f[length(f)-2],f[length(f)-1],f[length(f)], sep = "")
  # If extension aff then single file read, if txt then a list of file
  if (f_ext == "aff") {
    nfiles <- 1
    files <- c(file_name)
  } else if (f_ext == "txt") {
    con <- file(paste(getwd(),"/",file_name, sep = ""), "r")
    files <- readLines(con = con, n = -1)
    nfiles <- length(files)
    close(con)
  } else {
    return("Error: File extension must be .aff or .txt")
  }
  # Format momentum
  momentum_asint <- as.integer(nucleon_momentum)
  nucleon_momentum <- format_momentum(nucleon_momentum)
  # Create the parity projector
  parity_projector <- (id + parity*g0)/2
  # Read source coordinates of all files
  source_coord <- array(data = NA, dim = c(nfiles,4))
  for (i in 1:nfiles) {
    source_coord[i,] <- get_source_coord(files[i])
  }
  nread <- LT*4*4
  corr <- array(NA,dim=c(nfiles,nread))
  if (read_from_disk) {
    for (i in 1:nfiles) {
      # Set the files to read from
      file_to_read_t1 <- paste(read_write_dir,files[i],".",nucleon_type,"-",nucleon_type,".Gi_",Gi,".Gf_",Gf,".t1.px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
      file_to_read_t2 <- paste(read_write_dir,files[i],".",nucleon_type,"-",nucleon_type,".Gi_",Gi,".Gf_",Gf,".t2.px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
      # Read
      a <- read.csv(file = file_to_read_t1, header = FALSE)[,1]
      b <- read.csv(file = file_to_read_t2, header = FALSE)[,1]
      # Calculate correlator
      corr[i,] <- (a+b)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    }
  }else {
    for (i in 1:nfiles) {
      fwd <- paste(file_dir,files[i],sep="")
      # Generate key
      key <- generate_key_NN(nucleon_type=nucleon_type, source_coord=source_coord[i,], Gi=Gi, Gf=Gf, nucleon_momentum=nucleon_momentum)
      # Read keys from aff file
      a <- aff_read_key(filename=fwd, key=key[1], key_length=nread)
      b <- aff_read_key(filename=fwd, key=key[2], key_length=nread)
      if (write_to_disk) {
        # Set files to write to
        file_to_write_t1 <- paste(read_write_dir,files[i],".",nucleon_type,"-",nucleon_type,".Gi_",Gi,".Gf_",Gf,".t1.px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
        file_to_write_t2 <- paste(read_write_dir,files[i],".",nucleon_type,"-",nucleon_type,".Gi_",Gi,".Gf_",Gf,".t2.px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
        # Write
        write.table(a, file = file_to_write_t1, row.names = FALSE, col.names = FALSE)
        write.table(b, file = file_to_write_t2, row.names = FALSE, col.names = FALSE)
      }
      # Calculate correlator
      corr[i,] <- (a+b)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    }
  }
  # Parity projection and spin trace
  traced_corr <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    for (t in 1:LT) {
      for (spi in 1:4) {
        for (spj in 1:4) {
          traced_corr[i,t] <- traced_corr[i,t] + parity_projector[spi,spj] * corr[i,((t-1)*4 + spj-1)*4+spi]
        }
      }
    }
  }
  # Anti-periodic boundary conditions to get the final time dependent correlator
  corr_t <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    ti <- as.integer(source_coord[i,1])
    for (t in 0:LT-1) {
      delt <- t - ti
      if (delt >= 0 & delt <= (LT-ti-1))  {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*delt/LT)*traced_corr[i,t+1]
      } else if (delt < 0 & delt >= (-ti)) {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*(LT+delt)/LT)*traced_corr[i,t+1]
      }
    }
  }
  return(corr_t)
}

# Function that converts correlator arrays to cf objects of the Hadron library
corr_to_cf <- function(corr) {
  dims <- dim(corr)
  len <- length(dims)
  if (len == 2) {
    nmeasure <- dims[1]
    LT <- dims[2]
    dc <- array(NA, dim=c(nmeasure,LT,c(1,1)))
    for (i in 1:nmeasure) {
      for (j in 1:LT) {
        dc[i,j,1,1] <- corr[i,j]
      }
    }
  }
  if (len == 3) {
    LT <- dims[3]
    nmeas <- dims[1]*dims[2]
    dc <- array(NA, dim=c(nmeas,LT,c(1,1)))
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:LT) {
          dc[(i-1)*dims[2]+j,k,1,1] <- corr[i,j,k]
        }
      }
    }
  }else if (len == 4) {
    LT <- dims[4]
    nmeas <- dims[1]*dims[2]*dims[3]
    dc <- array(NA, dim=c(nmeas,LT,c(1,1)))
    for (i in 1:dims[1]) {
      for (j in 1:dims[2]) {
        for (k in 1:dims[3]) {
          for (l in 1:LT) {
            dc[((i-1)*dims[2]+j-1)*dims[3]+k,l,1,1] <- corr[i,j,k,l]
          }
        }
      }
    }
  }
  
  cmat <- raw_cf_meta(cf = raw_cf(), nrObs = 1, Time = LT, nrStypes = 1, dim = c(1,1))
  cmat <- raw_cf_data(cmat, data = dc)
  cmat <- raw_cf_to_cf(cmat, component = c(1,1))
  return(cmat)
}
#   if (len == 3) {
#     LT <- dims[3]
#     nmeas <- dims[1]*dims[2]
#     dc <- array(NA, dim=c(nmeas,LT,c(1,1)))
#     for (i in 1:dims[1]) {
#       for (j in 1:dims[2]) {
#         for (k in 1:LT) {
#           dc[(i-1)*dims[2]+j,k,1,1] <- corr[i,j,k]
#         }
#       }
#     }
#   }else if (len == 4) {
#     LT <- dims[4]
#     nmeas <- dims[1]*dims[2]*dims[3]
#     dc <- array(NA, dim=c(nmeas,LT,c(1,1)))
#     for (i in 1:dims[1]) {
#       for (j in 1:dims[2]) {
#         for (k in 1:dims[3]) {
#           for (l in 1:LT) {
#             dc[((i-1)*dims[2]+j-1)*dims[3]+k,l,1,1] <- corr[i,j,k,l]
#           }
#         }
#       }
#     }
#   }
#   
#   cmat <- raw_cf_meta(cf = raw_cf(), nrObs = 1, Time = LT, nrStypes = 1, dim = c(1,1))
#   cmat <- raw_cf_data(cmat, data = dc)
#   cmat <- raw_cf_to_cf(cmat, component = c(1,1))
#   return(cmat)
# }

# Function to generate file name to read from or write to disk
generate_file_name <- function(file_name, diagram, diagram_type, Gc, Gi, Gf, nucleon_momentum, insertion_momentum, nsample, sample) {
  if (diagram_type == "W") {
    io_file_name <- paste(file_name,".",diagram[1],".QX",insertion_momentum[1],"_QY",insertion_momentum[2],"_QZ",insertion_momentum[3],".sample",sample,".Gc_",Gc,".Gi_",
                          Gi,".Gf_",Gf,".t",diagram[2],".px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
  }else {
    io_file_name <- paste(file_name,".",diagram[1],".QX",insertion_momentum[1],"_QY",insertion_momentum[2],"_QZ",insertion_momentum[3],".nsample",nsample,".Gc_",Gc,".Gi_",
                          Gi,".Gf_",Gf,".t",diagram[2],".px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],".txt",sep="")
  }
  return(io_file_name)
}

# Calculate correlators for NjjN
NjjN_corr <- function(file_name, file_dir = "/home/aniket/projects/npv_data/cA211a.30.32/", read_write_dir = NA, 
                      nucleon_type, insertion_type, set_dtype = 0, set_ttype = 0, diagram_type = 0, Gc = 0, Gi = "Cg5", Gf = "Cg5", 
                      nucleon_momentum = c(0,0,0), insertion_momentum = c(0,0,0), nsample = 1, samples = 1, 
                      sample_mean = TRUE, LT = 64, parity = +1, read_from_disk = FALSE, write_to_disk = FALSE) {
  # Split file name into characters and determine the extension
  f <- strsplit(file_name, split = "")[[1]]
  f_ext <- paste(f[length(f)-2],f[length(f)-1],f[length(f)], sep = "")
  # If the extension is txt then read all the files to read, if aff then a single file
  if (f_ext == "aff") {
    nfiles <- 1
    files <- c(file_name)
  } else if (f_ext == "txt") {
    con <- file(paste(getwd(),"/",file_name, sep = ""), "r")
    files <- readLines(con = con, n = -1)
    close(con)
    nfiles <- length(files)
  } else {
    return("Error: File extension must be .aff or .txt")
  }
  # Define the parity projection operator
  parity_projector <- (id + parity*g0)/2
  # Get source coordinates from the files
  source_coord <- array(data = NA, dim = c(nfiles,4))
  for (i in 1:nfiles) {
    source_coord[i,] <- get_source_coord(files[i])
  }
  momentum_int <- as.numeric(nucleon_momentum)
  nucleon_momentum <- format_momentum(nucleon_momentum)
  nread <- LT*4*4
  if (set_dtype == 0 | set_ttype == 0){
    diagrams <- gen_diag(nucleon = nucleon_type, insertion = insertion_type, dtype = diagram_type, set_dtype = 0, set_ttype = 0)
  }else {
    diagrams <- unique(gen_diag(nucleon = nucleon_type, insertion = insertion_type, dtype = diagram_type, set_dtype = set_dtype, set_ttype = set_ttype))
  }
  if (sample_mean) {
  corr <- array(0,dim=c(nfiles,nread))
  #diagrams <- gen_diag(nucleon = nucleon_type, insertion = insertion_type, dtype = diagram_type)
  for (i in 1:nfiles) {
    if (read_from_disk) {
      for (d in diagrams) {
        sample_read <- array(0,dim=c(samples,nread))
        sample_mean <- array(0,dim=nread)
        for (s in 1:samples) {
          file_to_read <- paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d, diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, 
                                                                  nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s),sep="")
          sample_read[s,] <- read.csv(file = file_to_read, header = FALSE)[,1]
        }
        for (t in 1:nread) sample_mean[t] <- mean(sample_read[,t])
        corr[i,] <- corr[i,] + sample_mean
      }
      corr[i,] <- corr[i,]*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
    }else {
      keys <- c()
      if (write_to_disk) files_to_write <- c()
      for (d in diagrams) {
        sample_read <- array(0,dim=c(samples,nread))
        sample_mean <- array(0,dim=nread)
        for (s in 1:samples) {
          keys <- c(keys,generate_key(ts = d, sc = source_coord[i,], mom = nucleon_momentum, jmom = as.integer(insertion_momentum), nsample = as.integer(nsample),
                                      sample = as.integer(s-1), Gi = Gi, Gf = Gf, dtype = diagram_type, vtype = Gc))
          if (write_to_disk) files_to_write <- c(files_to_write,paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d,diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s),sep=""))
        }
      }
      print(keys)
      file_to_read <- paste(file_dir,files[i],sep="")
      tempread <- aff_read_key_list(filename = file_to_read, key_list = keys, key_length = nread)
      for (l in 1:length(keys)) {
        corr[i,] <- corr[i,] + tempread[((l-1)*nread+1):(l*nread)]
        if (write_to_disk) write.table(tempread[((l-1)*nread+1):(l*nread)], file = files_to_write[l], row.names = FALSE, col.names = FALSE)
      }
      corr[i,] <- (corr[i,]/samples)*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
    }
  }
  #   for (d in diagrams) {
  #     sample_read <- array(0,dim=c(samples,nread))
  #     sample_mean <- array(0,dim=nread)
  #     if (read_from_disk) {
  #       for (s in 1:samples) {
  #         file_to_read <- paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d, diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, 
  #                                                                 nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s),sep="")
  #         sample_read[s,] <- read.csv(file = file_to_read, header = FALSE)[,1]
  #       }
  #     }
  #     for (s in 1:samples) {
  #       if (read_from_disk) {
  #         file_to_read <- paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d, diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, 
  #                                            nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s),sep="")
  #         sample_read[s,] <- read.csv(file = file_to_read, header = FALSE)[,1]
  #       }else {
  #         key <- generate_key(ts = d, sc = source_coord[i,], mom = nucleon_momentum, jmom = as.integer(insertion_momentum), nsample = as.integer(nsample),
  #                             sample = as.integer(s-1), Gi = Gi, Gf = Gf, dtype = diagram_type, vtype = Gc)
  #         file_to_read <- paste(file_dir,files[i],sep="")
  #         sample_read[s,] <- aff_read_key(filename = file_to_read, key = key, key_length = nread)
  #         if (write_to_disk) {
  #           file_to_write <- paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d, diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, 
  #                                                                            nucleon_momentum = nucleon_momentum, insertion_momentum = insertion_momentum, nsample = nsample, sample = s),sep="")
  #           write.table(sample_read[s,], file = file_to_write, row.names = FALSE, col.names = FALSE)
  #         }
  #       }
  #     }
  #     for (t in 1:nread) sample_mean[t] <- mean(sample_read[,t])
  #     corr[i,] <- corr[i,] + sample_mean
  #   }
  #   corr[i,] <- corr[i,]*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
  # }
  #   fwd <- paste(fdir,files[i],sep="")
  #   for (j in 1:nmoms) {
  #     for (k in 1:jmoms) {
  #       sread <- array(0,dim=c(samples,nread))
  #       ssum <- array(0,dim=nread)
  #       for (s in 1:samples) {
  #       keys <- generate_key(ts = ts, source_coord = source_coord[i,], mom = nmomlist[j,], jmom = jmomlist[k,], nsample = as.integer(nsample), 
  #                            sample = as.integer(s-1), Gi = Gi, Gf = Gf, dtype = dtype, vtype = Gc)
  #       lkey <- length(keys)
  #       tempread <- aff_read_key_list(filename = fwd, key_list = keys, key_length = nread)
  #       for (l in 1:lkey) {
  #         sread[s,] <- sread[s,] + tempread[((l-1)*nread+1):(l*nread)]
  #       }
  #       }
  #       for (t in 1:nread) ssum[t] <- mean(sread[,t])
  #       mm <- as.numeric(nmomlist[j,])
  #       corr[i,j,k,] <- ssum*exp(-1.i*(mm[1]*source_coord[i,2]+mm[2]*source_coord[i,3]+mm[3]*source_coord[i,4]))
  #     }  
  #   }
  # }
  
  # Parity projection and spin trace
  traced_corr <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    for (t in 1:LT) {
      for (spi in 1:4) {
        for (spj in 1:4) {
          traced_corr[i,t] <- traced_corr[i,t] + parity_projector[spi,spj] * corr[i,((t-1)*4 + spj-1)*4+spi]
        }
      }
    }
  }
    
  # Anti-periodic boundary condition to get the final NjjN correlator
  corr_t <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    ti <- as.integer(source_coord[i,1])
    for (t in 0:LT-1) {
      delt <- t - ti
      if (delt >= 0 & delt <= (LT-ti-1))  {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*delt/LT)*traced_corr[i,t+1]
      } else if (delt < 0 & delt >= (-ti)) {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*(LT+delt)/LT)*traced_corr[i,t+1]
      }
    }
  }
  }else {
    corr <- array(0,dim=c(nfiles,samples,nread))
    #diagrams <- gen_diag(nucleon = nucleon_type, insertion = insertion_type, dtype = diagram_type)
    for (i in 1:nfiles) {
      if (read_from_disk) {
        for (d in diagrams) {
          for (s in 1:samples) {
            file_to_read <- paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d, diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, 
                                                                    nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s),sep="")
            corr[i,s,] <- corr[i,s,] + read.csv(file = file_to_read, header = FALSE)[,1]
          }
        }
        corr[i,,] <- corr[i,,]*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
      }else {
        keys <- c()
        if (write_to_disk) files_to_write <- c()
        for (d in diagrams) {
          sample_read <- array(0,dim=c(samples,nread))
          sample_mean <- array(0,dim=nread)
          for (s in 1:samples) {
            keys <- c(keys,generate_key(ts = d, sc = source_coord[i,], mom = nucleon_momentum, jmom = as.integer(insertion_momentum), nsample = as.integer(nsample),
                                        sample = as.integer(s-1), Gi = Gi, Gf = Gf, dtype = diagram_type, vtype = Gc))
            if (write_to_disk) files_to_write <- c(files_to_write,paste(read_write_dir,generate_file_name(file_name = files[i], diagram = d,diagram_type = diagram_type, Gc = Gc, Gi = Gi, Gf = Gf, nucleon_momentum = nucleon_momentum, insertion_momentum = as.character(insertion_momentum), nsample = nsample, sample = s-1),sep=""))
          }
        }
        file_to_read <- paste(file_dir,files[i],sep="")
        tempread <- aff_read_key_list(filename = file_to_read, key_list = keys, key_length = nread)
        for (l in 1:length(keys)) {
          s <- l%%samples+1
          corr[i,s,] <- corr[i,s,] + tempread[((l-1)*nread+1):(l*nread)]
          if (write_to_disk) write.table(tempread[((l-1)*nread+1):(l*nread)], file = files_to_write[l], row.names = FALSE, col.names = FALSE)
        }
        corr[i,,] <- corr[i,,]*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
      }
    }
    traced_corr <- array(0, dim = c(nfiles,samples,LT))
    for (i in 1:nfiles) {
      for (s in 1:samples){
      for (t in 1:LT) {
        for (spi in 1:4) {
          for (spj in 1:4) {
            traced_corr[i,s,t] <- traced_corr[i,s,t] + parity_projector[spi,spj] * corr[i,s,((t-1)*4 + spj-1)*4+spi]
          }
        }
      }
      }
    }
    
    # Anti-periodic boundary condition to get the final NjjN correlator
    corr_t <- array(0, dim = c(nfiles,samples,LT))
    for (i in 1:nfiles) {
      ti <- as.integer(source_coord[i,1])
      for (s in 1:samples) {
      for (t in 0:LT-1) {
        delt <- t - ti
        if (delt >= 0 & delt <= (LT-ti-1))  {
          corr_t[i,s,(LT+delt)%%LT+1] <- exp(1.i*3*pi*delt/LT)*traced_corr[i,s,t+1]
        } else if (delt < 0 & delt >= (-ti)) {
          corr_t[i,s,(LT+delt)%%LT+1] <- exp(1.i*3*pi*(LT+delt)/LT)*traced_corr[i,s,t+1]
        }
      }
      }
    }
  }
  return(corr_t)
}

# Function to bootstrap correlation function and save the mean and error
bootstrap_cf_write <- function(cmat, R = 400, l = 2, seed = 1234, sim = "geom", endcorr = TRUE, file_to_write) {
  cmat <- bootstrap.cf(cmat, boot.R = R, boot.l = l, seed = seed, sim = sim, endcorr = endcorr)
  dat <- cmat$cf0
  se <- cmat$tsboot.se
  idat <- cmat$icf0
  ise <- cmat$itsboot.se
  df <- data.frame(dat,se,idat,ise)
  write.table(df,file_to_write)
}

read_and_bootstrap <- function(file_to_read,file_to_write) {
  corr <- readRDS(file_to_read)
  cmat <- corr_to_cf(corr = corr)
  bootstrap_cf_write(cmat = cmat, file_to_write = file_to_write)
}

plot.new <- function(x, neg.vec = rep(1, times = length(cf$cf0)), x.fac = 0., y.fac = 0., rep = FALSE, ...) {
  cf <- x
  stopifnot(inherits(cf, 'cf_boot'))
  stopifnot(inherits(cf, 'cf_meta'))
  
  val <- cf$cf0
  err <- cf$tsboot.se
  
  if(!cf$symmetrised){
    tmax <- cf$Time - 1
  } else {
    tmax <- cf$Time / 2
  }
  
  df <- data.frame(t = rep(c(0:tmax), times = length(val)/(tmax+1)),
                   CF = val,
                   Err = err)
  
  plotwitherror(x = df$t + x.fac, y = neg.vec * df$CF + y.fac, dy = df$Err, rep = rep, ...)
  
  return(invisible(df))
}

plot_mult_new <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", xfac = 0., leg, legpos = NA, cex = 2.5, cex.lab=1.0, cex.axis=1.0, axes=TRUE, cex.leg = cex, showleg = TRUE, lwd = 1) {
  plot_files <- list(...)
  cols <- c("dodgerblue3","firebrick1","chartreuse3","brown","orange")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot.new(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],x.fac=xfac[1],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot.new(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,x.fac=xfac[i],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  grid(col="gray38", lwd=1.)
  if (!is.na(title)) title(title)
  if (showleg == TRUE) {
    if (is.na(legpos)) {
      if (is.na(ylim)) legpos <- c(17,y0)
      else legpos <- c(17,ylim[2])
    }
    legend(legpos[1],legpos[2],leg,col=cols_used,pch=pch_used,cex=cex.leg,bg="white",box.col="gray38")
  }
}

plot_mult_new_nogrid <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", xfac = 0., leg, legpos = NA, cex = 2.5, cex.lab=1.0, cex.axis=1.0, axes=TRUE, cex.leg = cex, showleg = TRUE, lwd = 1) {
  plot_files <- list(...)
  cols <- c("dodgerblue3","firebrick1","chartreuse3","brown","orange")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot.new(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],x.fac=xfac[1],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot.new(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,x.fac=xfac[i],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (showleg == TRUE) {
    if (is.na(legpos)) {
      if (is.na(ylim)) legpos <- c(17,y0)
      else legpos <- c(17,ylim[2])
    }
    legend(legpos[1],legpos[2],leg,col="transparent",pch=pch_used,cex=cex.leg,bg="transparent",box.col="transparent") #"gray38")
  }
}

plot_mult_new_nogrid_noabline <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", xfac = 0., leg, legpos = NA, cex = 2.5, cex.lab=1.0, cex.axis=1.0, axes=TRUE, cex.leg = cex, showleg = TRUE, lwd = 1) {
  plot_files <- list(...)
  cols <- c("dodgerblue3","firebrick1","chartreuse3","brown","orange")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot.new(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],x.fac=xfac[1],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot.new(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,x.fac=xfac[i],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  if (!is.na(title)) title(title)
  if (showleg == TRUE) {
    if (is.na(legpos)) {
      if (is.na(ylim)) legpos <- c(17,y0)
      else legpos <- c(17,ylim[2])
    }
    legend(legpos[1],legpos[2],leg,col="transparent",pch=pch_used,cex=cex.leg,bg="transparent",box.col="transparent") #"gray38")
  }
}

plot_mult_new_nogrid_transparent_leg <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", xfac = 0., leg, legpos = NA, cex = 2.5, cex.lab=1.0, cex.axis=1.0, axes=TRUE, cex.leg = cex, showleg = TRUE, lwd = 1) {
  plot_files <- list(...)
  cols <- c("dodgerblue3","firebrick1","chartreuse3","brown","orange")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot.new(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],x.fac=xfac[1],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot.new(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,x.fac=xfac[i],cex=cex,cex.lab=cex.lab,cex.axis=cex.axis,axes=axes,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (showleg == TRUE) {
    if (is.na(legpos)) {
      if (is.na(ylim)) legpos <- c(17,y0)
      else legpos <- c(17,ylim[2])
    }
    legend(legpos[1],legpos[2],leg,col="transparent",pch=pch_used,cex=cex.leg,bg=alpha( col="white", a=0.75 ),box.col=alpha( col="white", a=0.75 )) #"gray38")
  }
}

plot_mult_new2 <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", leg, legpos = NA, cex = 2.5, cex.leg = cex, lwd = 1) {
  plot_files <- list(...)
  cols <- c("navy","royalblue3","deepskyblue2","maroon","red","tomato","olivedrab","green3","yellowgreen")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,cex=cex,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (is.na(legpos)) {
    if (is.na(ylim)) legpos <- c(17,y0)
    else legpos <- c(17,ylim[2])
  }
  grid(col="gray38", lwd=1.)
  legend(legpos[1],legpos[2],leg,col=cols_used,pch=pch_used,cex=cex.leg,bg="white",box.col="gray38")
}

plot_mult <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", leg, legpos = NA, cex = 2.5, y.intersp.leg=1.0, cex.leg = cex, cex.axis = cex, cex.lab = cex, lwd = 1) {
  plot_files <- list(...)
  cols <- c("darkorange","chartreuse3","deepskyblue1","firebrick1","darkgreen","navyblue","orange")
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab)
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,cex=cex,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (is.na(legpos)) {
    if (is.na(ylim)) legpos <- c(17,y0)
    else legpos <- c(17,ylim[2])
  }
  grid(col="gray38", lwd=1.)
  legend(legpos[1],legpos[2],leg,col=cols_used,pch=pch_used,cex=cex.leg,bg="white",box.col="gray38",y.intersp=y.intersp.leg)
}

plot_mult_nogrid <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", leg, legpos = NA, cex = 2.5, y.intersp.leg=1.0, cex.leg = cex, cex.axis = cex, cex.lab = cex, lwd = 1,
                             no_xlabels = FALSE, version = 1) {
  plot_files <- list(...)
  if (version == 1) {
    cols <- c("darkorange","chartreuse3","deepskyblue1","firebrick1","darkgreen","navyblue","orange")
  } else if (version == 2) {
    cols <- c("darkorange","chartreuse3","chartreuse4","deepskyblue1","firebrick1","darkgreen","navyblue","orange")
  } else {
    print("Error: Select either version=1 or version=2.")
    stop()
  }
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  if (no_xlabels==FALSE) {
    plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab)
  } else {
    plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab, xaxt="n")
    axis(side=1, cex.axis=1.4, at=seq(from=0,to=20,by=5), labels=FALSE)
  }
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,cex=cex,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (is.na(legpos)) {
    if (is.na(ylim)) legpos <- c(17,y0)
    else legpos <- c(17,ylim[2])
  }
  legend(legpos[1],legpos[2],leg,col=cols_used,pch=pch_used,cex=cex.leg,bg="white",box.col="gray38",y.intersp=y.intersp.leg)
}

plot_single_nogrid <- function(..., title = NA, ylim = NA, xlim = c(0,20), xlab = "x", ylab = "y", leg, legpos = NA, cex = 2.5, y.intersp.leg=1.0, cex.leg = cex, cex.axis = cex, cex.lab = cex, lwd = 1,
                             no_xlabels = FALSE, version = 1) {
  plot_files <- list(...)
  if (version == 1) {
    cols <- c("firebrick1","firebrick1","firebrick1","firebrick1","firebrick1","firebrick1","firebrick1")
  } else if (version == 2) {
    cols <- c("firebrick1","firebrick1","firebrick1","firebrick1","firebrick1","firebrick1","firebrick1","firebrick1")
  } else {
    print("Error: Select either version=1 or version=2.")
    stop()
  }
  # pch_used <- seq(1,length(plot_files))
  pch_used <- rep(16,length(plot_files))
  if (is.na(ylim)) {
    mx <- c()
    mn <- c()
    for (i in c(1:length(plot_files))) {
      mx <- c(mx,max(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
      mn <- c(mn,min(plot_files[i][[1]]$cf0[(xlim[1]+1):(xlim[2]+1)]))
    }
    y0 <- max(mx)*1.2
    y1 <- min(mn)*0.8
    ylim <- c(y1,y0)
  }
  if (no_xlabels==FALSE) {
    plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab)
  } else {
    plot(plot_files[1][[1]],col=cols[1],ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,pch=pch_used[1],cex=cex,lwd=lwd, cex.axis=cex.axis, cex.lab=cex.lab, xaxt="n")
    axis(side=1, cex.axis=1.4, at=seq(from=0,to=20,by=5), labels=FALSE)
  }
  cols_used <- cols[1]
  for (i in 2:length(plot_files)) {
    plot(plot_files[i][[1]],col=cols[i],pch=pch_used[i],rep=TRUE,cex=cex,lwd=lwd)
    cols_used <- c(cols_used,cols[i])
  }
  abline(h=0)
  if (!is.na(title)) title(title)
  if (is.na(legpos)) {
    if (is.na(ylim)) legpos <- c(17,y0)
    else legpos <- c(17,ylim[2])
  }
  legend(legpos[1],legpos[2],leg,col=cols_used,pch=pch_used,cex=cex.leg,bg="white",box.col="gray38",y.intersp=y.intersp.leg)
}

plt.errband <- function (x, y, se, lwd = 1.5, col="gray", a=0.25) {
  polyx <- c( x, rev(x) )
  polyval <- c( (y + se), rev(y - se) )
  ## for a = 0.25 the opacity is 25%
  polycol <- alpha( col, a )
  
  polygon( x=polyx, y=polyval, col=polycol, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}

plt.errband1 <- function (x, y, se1, se2, lwd = 1.5, col="gray38", a=0.25) {
  polyx <- c( x, rev(x) )
  polyval <- c( (y + se2), rev(y - se1) )
  ## for a = 0.25 the opacity is 25%
  polycol <- alpha( col, a )
  
  polygon( x=polyx, y=polyval, col=polycol, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}

plt.errband2 <- function (x, y, se1, se2, lwd = 1.5, col="gray", a=0.25) {
  polyx <- c( (x + se2), rev(x - se1) )
  polyval <- c( y, rev(y) )
  ## for a = 0.25 the opacity is 25%
  polycol <- alpha( col, a )
  
  polygon( x=polyx, y=polyval, col=polycol, border=FALSE )
  lines( x, y, col=col, lwd=lwd )
}

calc_dm <- function(cfl,cf,tau) {
  dm <- cfl/cf - shift.cf(cfl,tau)/shift.cf(cf,tau)
  dm <- mul.cf(dm,(1/tau))
  return(dm)
}

calc_dm2 <- function(cfl,cf,tau) {
  cfl1 <- shift.cf(cfl,tau)
  cf1 <- shift.cf(cf,tau)
  r1 <- cfl$cf0/cf$cf0
  dr1 <- sqrt((cfl$tsboot.se/cf$cf0)^2+(cfl$cf0*cf$tsboot.se/cf$cf0^2)^2)
  r2 <- cfl1$cf0/cf1$cf0
  dr2 <- sqrt((cfl1$tsboot.se/cf1$cf0)^2+(cfl1$cf0*cf1$tsboot.se/cf1$cf0^2)^2)
  # r <- (r1 - r2)/tau
  # c <- cov((cfl/cf)$cf)
  # dr <- sqrt(dr1^2+dr2^2)/tau
  # return(data.frame(r,dr))
}

# Function to read individual diagrams given a filelist
read_diagrams <- function(file_name, file_dir, write_dir, diagrams) {
  
}

read_rds_bootstrap <- function(file_to_read) {
  a <- readRDS(file_to_read)
  b <- corr_to_cf(a)
  b <- bootstrap.cf(b)
  return(b)
}

srcavg <- function(corr,nsrc,totsrc) {
  n <- dim(corr)[1]
  LT <- dim(corr)[2]
  n <- n/totsrc
  corrsrc <- array(NA, dim=c(n,LT))
  for (i in 1:n) {
    for (t in 1:LT) {
      corrsrc[i,t] <- mean(corr[((i-1)*totsrc+1):((i-1)*totsrc+nsrc),t])
    }
  }
  cmatsrc <- corr_to_cf(corrsrc)
  cmatsrc <- bootstrap.cf(cmatsrc)
  return(cmatsrc)
}

srcavg.raw <- function(corr,nsrc,totsrc) {
  n <- dim(corr)[1]
  LT <- dim(corr)[2]
  n <- n/totsrc
  corrsrc <- array(NA, dim=c(n,LT))
  for (i in 1:n) {
    for (t in 1:LT) {
      corrsrc[i,t] <- mean(corr[((i-1)*totsrc+1):((i-1)*totsrc+nsrc),t])
    }
  }
  return(corrsrc)
}

oetavg <- function(corr,n_oet,tot_oet) {
  m <- dim(corr)[1]
  n <- dim(corr)[2]
  LT <- dim(corr)[3]
  n <- n/tot_oet
  corrsrc <- array(NA, dim=c(m,LT))
  for (i in 1:n) {
    for (j in 1:m) {
      for (t in 1:LT) {
        corrsrc[j,t] <- mean(corr[j,((i-1)*tot_oet+1):((i-1)*tot_oet+n_oet),t])
      }
    }
  }
  return(corrsrc)
  # cmatsrc <- corr_to_cf(corrsrc)
  # cmatsrc <- bootstrap.cf(cmatsrc)
  # return(cmatsrc)
}

sample_mean <- function(corr) {
  n <- dim(corr)
  corr1 <- array(NA,dim=c(n[1],n[3]))
  for (i in 1:n[1]) {
    for (j in 1:n[3]) {
      corr1[i,j] <- mean(corr[i,,j])
    }
  }
  return(corr1)
}

archiver <- function(file_name,archive_dir) {
  con <- file(paste(getwd(),"/",file_name, sep = ""), "r")
  files <- readLines(con = con, n = -1)
  close(con)
  nfiles <- length(files)
  src_files <- ""
  for (i in 1:nfiles) {
    config <- paste(strsplit(files[i],split="")[[1]][10:13],collapse="")
    src_files <- paste(src_files," ",files[i])
    if (paste(strsplit(files[i+1],split="")[[1]][10:13],collapse="") != config) {
      cmd <- paste("tar -cvf ",archive_dir,"njjn_fht.",config,".tar ", src_files, sep="" )
      system(cmd)
      src_files <- ""
    }
  }
}

# Generate Key for NjjN Correlator with strange quark operator insertion
generate_key_NjjN_strange <- function(nucleon_type, source_coord, Gi = "Cg5", Gf = "Cg5", nucleon_momentum, jmom, nsample, insertion_type, vtype) {
  if (nucleon_type == 'p' && insertion_type == 'sdds') {
    nkeys <- 1
    print(nkeys)
  } else if (nucleon_type == 'n' && insertion_type == 'suus') {
    nkeys <- 1
    print(nkeys)
  } else if (nucleon_type == 'p' && insertion_type == 'suus') {
    nkeys <- 2
    print(nkeys)
  } else if (nucleon_type == 'n' && insertion_type == 'sdds') {
    nkeys <- 2
    print(nkeys)
  }
  if ( nkeys == 1 ) {
    # Set the propagator ordering based on nucleon type
    if (nucleon_type == 'p'){ 
      k <- 'fu-bddd-fu'
    }  else if (nucleon_type == 'n') { 
      k <- 'fd-duuu-fd'
    }  else {
      print("nucleon_type must be 'n' or 'p'")
      return(NA)
    }
    # Write the keys for t1 and t2
    key1a <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t1/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key1b <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key <- c( key1a, key1b )
  } else {
    # Set the propagator ordering based on nucleon type
    if (nucleon_type == 'p'){ 
      k1 <- 'buuu-fd-fu'
      k2 <- 'fu-fd-buuu'
    }  else if (nucleon_type == 'n') { 
      k1 <- 'bddd-fu-fd'
      k2 <- 'fd-fu-bddd'
    }  else {
      print("nucleon_type must be 'n' or 'p'")
      return(NA)
    }
    # Write the keys for t1 and t2
    key1a <- paste("/N-qbGqqbGq-N/",k1,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t1/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key1b <- paste("/N-qbGqqbGq-N/",k1,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key2a <- paste("/N-qbGqqbGq-N/",k2,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t1/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key2b <- paste("/N-qbGqqbGq-N/",k2,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/nsample",nsample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    key <- c( key1a, key1b, key2a, key2b )
  }
  return(key)
}

# Calculate NjjN correlator for remaining strange quark insertion cases
# The function NjjN_corr_niko() was renamed and now is called NjjN_corr_strange()
NjjN_corr_strange <- function(file_name, file_dir = "/home/aniket/projects/npv_data/cA211a.30.32/", nucleon_type, Gc = 0, Gi = "Cg5", Gf = "Cg5", 
                           nucleon_momentum = c(0,0,0), LT = 64, parity = +1, read_from_disk = FALSE, write_to_disk = FALSE, read_write_dir = NA,
                           insertion_momentum = c(0,0,0), insertion_type, nsample = 1) {
  # Split the file name to get the extension
  f <- strsplit(file_name, split = "")[[1]]
  f_ext <- paste(f[length(f)-2],f[length(f)-1],f[length(f)], sep = "")
  # If extension aff then single file read, if txt then a list of file
  if (f_ext == "aff") {
    nfiles <- 1
    files <- c(file_name)
  } else if (f_ext == "txt") {
    con <- file(paste(getwd(),"/",file_name, sep = ""), "r")
    files <- readLines(con = con, n = -1)
    nfiles <- length(files)
    close(con)
  } else {
    return("Error: File extension must be .aff or .txt")
  }
  # Format momentum
  momentum_asint <- as.integer(nucleon_momentum)
  nucleon_momentum <- format_momentum(nucleon_momentum)
  # Create the parity projector
  parity_projector <- (id + parity*g0)/2
  # Read source coordinates of all files
  source_coord <- array(data = NA, dim = c(nfiles,4))
  for (i in 1:nfiles) {
    source_coord[i,] <- get_source_coord(files[i])
  }
  nread <- LT*4*4
  corr <- array(NA,dim=c(nfiles,nread))
  for (i in 1:nfiles) {
    fwd <- paste(file_dir,files[i],sep="")
    # Generate key
    key <- generate_key_NjjN_strange(nucleon_type=nucleon_type, source_coord=source_coord[i,], Gi=Gi, Gf=Gf, nucleon_momentum=nucleon_momentum, jmom = as.integer(insertion_momentum), nsample = as.integer(nsample), insertion_type = insertion_type, vtype=Gc)
    print(key)
    if (nucleon_type == 'p' && insertion_type == 'sdds') {
      # Read keys from aff file
      a <- aff_read_key(filename=fwd, key=key[1], key_length=nread)
      b <- aff_read_key(filename=fwd, key=key[2], key_length=nread)
      # Calculate correlator
      corr[i,] <- (a+b)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    } else if (nucleon_type == 'n' && insertion_type == 'suus') {
      # Read keys from aff file
      a <- aff_read_key(filename=fwd, key=key[1], key_length=nread)
      b <- aff_read_key(filename=fwd, key=key[2], key_length=nread)
      # Calculate correlator
      corr[i,] <- (a+b)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    } else if (nucleon_type == 'p' && insertion_type == 'suus') {
      # Read keys from aff file
      a <- aff_read_key(filename=fwd, key=key[1], key_length=nread)
      b <- aff_read_key(filename=fwd, key=key[2], key_length=nread)
      c <- aff_read_key(filename=fwd, key=key[3], key_length=nread)
      d <- aff_read_key(filename=fwd, key=key[4], key_length=nread)
      # Calculate correlator
      corr[i,] <- (a+b+c+d)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    } else if (nucleon_type == 'n' && insertion_type == 'sdds') {
      # Read keys from aff file
      a <- aff_read_key(filename=fwd, key=key[1], key_length=nread)
      b <- aff_read_key(filename=fwd, key=key[2], key_length=nread)
      c <- aff_read_key(filename=fwd, key=key[3], key_length=nread)
      d <- aff_read_key(filename=fwd, key=key[4], key_length=nread)
      # Calculate correlator
      corr[i,] <- (a+b+c+d)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
    }
  }
  # Parity projection and spin trace
  traced_corr <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    for (t in 1:LT) {
      for (spi in 1:4) {
        for (spj in 1:4) {
          traced_corr[i,t] <- traced_corr[i,t] + parity_projector[spi,spj] * corr[i,((t-1)*4 + spj-1)*4+spi]
        }
      }
    }
  }
  # Anti-periodic boundary conditions to get the final time dependent correlator
  corr_t <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    ti <- as.integer(source_coord[i,1])
    for (t in 0:LT-1) {
      delt <- t - ti
      if (delt >= 0 & delt <= (LT-ti-1))  {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*delt/LT)*traced_corr[i,t+1]
      } else if (delt < 0 & delt >= (-ti)) {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*(LT+delt)/LT)*traced_corr[i,t+1]
      }
    }
  }
  return(corr_t)
}

generate_key_NjjN_light_W <- function(nucleon_type, source_coord, Gi = "Cg5", Gf = "Cg5", nucleon_momentum, jmom, sample, insertion_type, vtype, colix, single_t) {
  # Set the propagator ordering based on nucleon type
  if (nucleon_type == 'p'){
    # For W type diagrams the sub-key for puup is 'wuu-fd-wuu'
    k <- 'wuu-fd-wuu'
  } else if (nucleon_type == 'n') {
    # For W type diagrams the sub-key for nddn is 'wdd-fu-wdd'
    k <- 'wdd-fu-wdd'
  } else {
    print("nucleon_type must be 'n' or 'p'")
    return(NA)
  }
  
  if ( single_t == TRUE ) {
    # create aff keys for crossed color indices (colix == TRUE) and non-color crossed indices (colix == FALSE)
    if ( colix == TRUE ) {
      # Write the keys for t1 and t2
      key1a <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t1/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
      key1b <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    } else {
      # Write the keys for t3 and t4
      key1a <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t3/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
      key1b <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t4/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    }
    key <- c( key1a, key1b )
  } else if ( single_t == FALSE ) {
    # create aff keys for crossed color indices (colix == TRUE) and non-color crossed indices (colix == FALSE)
    if ( colix == TRUE ) {
      # Write the keys for ( t1 + t2 )
      key <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t1_2/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    } else {
      # Write the keys for ( t3 + t4 )
      key <- paste("/N-qbGqqbGq-N/",k,"/T",source_coord[1],"_X",source_coord[2],"_Y",source_coord[3],"_Z",source_coord[4],"/QX",jmom[1],"_QY",jmom[2],"_QZ",jmom[3],"/sample",sample,"/Gc_",vtype,"/Gf_",Gf,"/Gi_",Gi,"/t3_4/px",nucleon_momentum[1],"py",nucleon_momentum[2],"pz",nucleon_momentum[3],sep="")
    }
  } else {
    print("'single_t' is a logical variable that is either 'TRUE' or 'FALSE'")
    return(NA)
  }
  
  return(key)
}

# if (sample_mean) {
#   corr <- array(0,dim=c(nfiles,nread))
#   #diagrams <- gen_diag(nucleon = nucleon_type, insertion = insertion_type, dtype = diagram_type)
#   for (i in 1:nfiles) {
#     keys <- c()
#     for (s in 1:samples) {
#       keys <- c(keys,generate_key(ts = d, sc = source_coord[i,], mom = nucleon_momentum, jmom = as.integer(insertion_momentum), nsample = as.integer(nsample),
#                                   sample = as.integer(s-1), Gi = Gi, Gf = Gf, dtype = diagram_type, vtype = Gc))
#     }
# 
#     print(keys)
#     file_to_read <- paste(file_dir,files[i],sep="")
#     tempread <- aff_read_key_list(filename = file_to_read, key_list = keys, key_length = nread)
#     for (l in 1:length(keys)) {
#       corr[i,] <- corr[i,] + tempread[((l-1)*nread+1):(l*nread)]
#     }
#     corr[i,] <- (corr[i,]/samples)*exp(-1.i*(momentum_int[1]*source_coord[i,2]+momentum_int[2]*source_coord[i,3]+momentum_int[3]*source_coord[i,4]))
#   }
# #===========================

# Calculate NjjN correlator for light quark insertion and W type diagram
NjjN_corr_light_W <- function(file_name, file_dir = "/home/aniket/projects/npv_data/cA211a.30.32/", nucleon_type, Gc = 0, Gi = "Cg5", Gf = "Cg5", 
                              nucleon_momentum = c(0,0,0), LT = 64, parity = +1, read_from_disk = FALSE, write_to_disk = FALSE, read_write_dir = NA,
                              insertion_momentum = c(0,0,0), insertion_type, samples = 1, colix, single_t) {
  # Split the file name to get the extension
  f <- strsplit(file_name, split = "")[[1]]
  f_ext <- paste(f[length(f)-2],f[length(f)-1],f[length(f)], sep = "")
  # If extension aff then single file read, if txt then a list of file
  if (f_ext == "aff") {
    nfiles <- 1
    files <- c(file_name)
  } else if (f_ext == "txt") {
    con <- file(paste(getwd(),"/",file_name, sep = ""), "r")
    files <- readLines(con = con, n = -1)
    nfiles <- length(files)
    close(con)
  } else {
    return("Error: File extension must be .aff or .txt")
  }
  # Format momentum
  momentum_asint <- as.integer(nucleon_momentum)
  nucleon_momentum <- format_momentum(nucleon_momentum)
  # Create the parity projector
  parity_projector <- (id + parity*g0)/2
  # Read source coordinates of all files
  source_coord <- array(data = NA, dim = c(nfiles,4))
  for (i in 1:nfiles) {
    source_coord[i,] <- get_source_coord(files[i])
  }
  
  nread <- LT*4*4
  corr <- array(0,dim=c(nfiles,nread))
  for (i in 1:nfiles) {
    keys <- c()
    for (s in 1:samples) {
      # Generate keys
      keys <- c(keys, generate_key_NjjN_light_W(nucleon_type=nucleon_type, source_coord=source_coord[i,], Gi=Gi, Gf=Gf, nucleon_momentum=nucleon_momentum, jmom = as.integer(insertion_momentum),
                                                sample = as.integer(s-1), insertion_type = insertion_type, vtype = Gc, colix = colix, single_t = single_t))
    }
    print(keys)
    
    file_to_read <- paste(file_dir,files[i],sep="")
    tmp_read <- aff_read_key_list(filename = file_to_read, key_list = keys, key_length = nread)
    for (l in 1:length(keys)) {
      corr[i,] <- corr[i,] + tmp_read[((l-1)*nread+1):(l*nread)]
    }
    corr[i,] <- (corr[i,]/samples)*exp(-1.i*(momentum_asint[1]*source_coord[i,2]+momentum_asint[2]*source_coord[i,3]+momentum_asint[3]*source_coord[i,4]))
  }
  
  # Parity projection and spin trace
  traced_corr <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    for (t in 1:LT) {
      for (spi in 1:4) {
        for (spj in 1:4) {
          traced_corr[i,t] <- traced_corr[i,t] + parity_projector[spi,spj] * corr[i,((t-1)*4 + spj-1)*4+spi]
        }
      }
    }
  }
  # Anti-periodic boundary conditions to get the final time dependent correlator
  corr_t <- array(0, dim = c(nfiles,LT))
  for (i in 1:nfiles) {
    ti <- as.integer(source_coord[i,1])
    for (t in 0:LT-1) {
      delt <- t - ti
      if (delt >= 0 & delt <= (LT-ti-1))  {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*delt/LT)*traced_corr[i,t+1]
      } else if (delt < 0 & delt >= (-ti)) {
        corr_t[i,(LT+delt)%%LT+1] <- exp(1.i*3*pi*(LT+delt)/LT)*traced_corr[i,t+1]
      }
    }
  }
  return(corr_t)
}

# function to remove im from NjjN_corr() data
remove_im <- function ( corr ) {
  corr1 <- as.character(corr)
  corr2 <- as.character(corr)
  for ( i in 1:length(corr) ) {
    a <- strsplit(corr1[i],"")[[1]]
    t <- which(a == ".")[2]-3
    corr1[i] <- substr(corr2[i], 1, t)
  }
  return( corr1 )
}

# function to remove re from NjjN_corr() data
remove_re <- function ( corr ) {
  corr1 <- as.character(corr)
  corr2 <- remove_im(corr)
  
  for ( i in 1:length(corr) ) {
    corr1[i] <- substr(corr1[i], (str_length(corr2[i]) + 1), (str_length(corr1[i]) - 1))
  }
  return( corr1 )
}

# read_file_dir <- function( nconf ) {
#   file_dir <- paste0("njjn_fht_strange_000",nconf,"/conf_000",nconf,"/")
#   return(file_dir)
# }
# read_file_list <- function( nconf ) {
#   file_list <- paste0("njjn_fht_strange_000",nconf,"/filelist.txt")
#   return(file_list)
# }

read_file_dir <- function( nconf ) {
  file_dir <- paste0("njjn_w_strange_000",nconf,"/conf_000",nconf,"/")
  return(file_dir)
}
read_file_list <- function( nconf ) {
  file_list <- paste0("njjn_w_strange_000",nconf,"/filelist.txt")
  return(file_list)
}
  
# function to calculate the NN correlator data for multiple gauge configurations using NN_corr()
# results are stored in .rds files
# print_NN_corr_list <- function ( ensemble_id="cA211.30.32", nucleon_type, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, T ) {
#   nconf <- nconf_f - nconf_i + 1
#   corrlist <- NA
#   for ( i in (nconf_i-1):(nconf_f-1) ) {
#     corr_tmp <- NN_corr(
#       file_name = read_file_list( i*4 ),
#       file_dir = read_file_dir( i*4 ),
#       nucleon_type = as.character( nucleon_type ),
#       parity = parity,
#       LT = T )
#     nsrc <- nrow(corr_tmp)
#     corrlist <- c(corrlist, colMeans(corr_tmp))
#   }
#   corrlist <- corrlist[-1]
#   t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
#   
#   if ( Re_Im == "re" ) {
#     corrlist <- cbind(t, remove_im(corrlist))
#   } else if ( Re_Im == "im" ) {
#     corrlist <- cbind(t, remove_re(corrlist))
#   } else if ( Re_Im == "re_im" ) {
#     corr_re <- remove_im(corrlist)
#     corr_im <- remove_re(corrlist)
#     corrlist <- cbind(t, corr_re, corr_im)
#   } else {
#     return("Error: Re_Im must be re, im or re_im")
#   }
#   colnames(corrlist) <- NULL
#   corrlist <- noquote(corrlist)
#   # print(corrlist)
#   saveRDS(corrlist, file = paste0(nucleon_type,nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".parity_",parity,".",Re_Im,".rds"))
# }

print_NN_corr_list <- function ( ensemble_id="cA211.30.32", nucleon_type, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, conf_list = c(), T ) {
  corrlist <- NA
  
  if ( is.null(conf_list) == TRUE ) {
    nconf <- nconf_f - nconf_i + 1
    for ( i in (nconf_i-1):(nconf_f-1) ) {
      corr_tmp <- NN_corr(
        file_name = read_file_list( i*4 ),
        file_dir = read_file_dir( i*4 ),
        nucleon_type = as.character( nucleon_type ),
        parity = parity,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  } else {
    nconf <- length(conf_list)
    for ( i in 1:nconf ) {
      corr_tmp <- NN_corr(
        file_name = read_file_list( conf_list[i] ),
        file_dir = read_file_dir( conf_list[i] ),
        nucleon_type = as.character( nucleon_type ),
        parity = parity,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  }
  corrlist <- corrlist[-1]
  t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
  
  if ( Re_Im == "re" ) {
    corrlist <- cbind(t, remove_im(corrlist))
  } else if ( Re_Im == "im" ) {
    corrlist <- cbind(t, remove_re(corrlist))
  } else if ( Re_Im == "re_im" ) {
    corr_re <- remove_im(corrlist)
    corr_im <- remove_re(corrlist)
    corrlist <- cbind(t, corr_re, corr_im)
  } else {
    return("Error: Re_Im must be re, im or re_im")
  }
  colnames(corrlist) <- NULL
  corrlist <- noquote(corrlist)
  # print(corrlist)
  saveRDS(corrlist, file = paste0(nucleon_type,nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".parity_",parity,".",Re_Im,".rds"))
}

# function to calculate the NjjN correlator data for multiple gauge configurations using NjjN_corr()
# results are stored in .rds files
# print_NjjN_corr_list <- function ( ensemble_id="cA211.30.32", nucleon_type , insertion_type, diagram_type, Gc, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, T ) {
#   nconf <- nconf_f - nconf_i + 1
#   corrlist <- NA
#   for ( i in (nconf_i-1):(nconf_f-1) ) {
#     corr_tmp <- NjjN_corr(
#       file_name = read_file_list( i*4 ),
#       file_dir = read_file_dir( i*4 ),
#       nucleon_type = as.character( nucleon_type ),
#       insertion_type = as.character( insertion_type ),
#       diagram_type = as.character( diagram_type ),
#       Gc = as.character( Gc ),
#       parity = parity,
#       LT = T )
#     nsrc <- nrow(corr_tmp)
#     corrlist <- c(corrlist, colMeans(corr_tmp))
#   }
#   corrlist <- corrlist[-1]
#   t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
#   
#   if ( Re_Im == "re" ) {
#     corrlist <- cbind(t, remove_im(corrlist))
#   } else if ( Re_Im == "im" ) {
#     corrlist <- cbind(t, remove_re(corrlist))
#   } else if ( Re_Im == "re_im" ) {
#     corr_re <- remove_im(corrlist)
#     corr_im <- remove_re(corrlist)
#     corrlist <- cbind(t, corr_re, corr_im)
#   } else {
#     return("Error: Re_Im must be re, im or re_im")
#   }
#   colnames(corrlist) <- NULL
#   corrlist <- noquote(corrlist)
#   # print(corrlist)
#   saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
# }

print_NjjN_corr_list <- function ( ensemble_id="cA211.30.32", nucleon_type , insertion_type, diagram_type, Gc, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, conf_list = c(), oet_samples = 1, T ) {
  corrlist <- NA
  
  if ( is.null(conf_list) == TRUE ) {
    nconf <- nconf_f - nconf_i + 1
    for ( i in (nconf_i-1):(nconf_f-1) ) {
      corr_tmp <- NjjN_corr(
        file_name = read_file_list( i*4 ),
        file_dir = read_file_dir( i*4 ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        diagram_type = as.character( diagram_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  } else {
    nconf <- length(conf_list)
    for ( i in 1:nconf ) {
      corr_tmp <- NjjN_corr(
        file_name = read_file_list( conf_list[i] ),
        file_dir = read_file_dir( conf_list[i] ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        diagram_type = as.character( diagram_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  }
  corrlist <- corrlist[-1]
  t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
  
  if ( Re_Im == "re" ) {
    corrlist <- cbind(t, remove_im(corrlist))
  } else if ( Re_Im == "im" ) {
    corrlist <- cbind(t, remove_re(corrlist))
  } else if ( Re_Im == "re_im" ) {
    corr_re <- remove_im(corrlist)
    corr_im <- remove_re(corrlist)
    corrlist <- cbind(t, corr_re, corr_im)
  } else {
    return("Error: Re_Im must be re, im or re_im")
  }
  colnames(corrlist) <- NULL
  corrlist <- noquote(corrlist)
  # print(corrlist)
  
  if (diagram_type == "W") {
    saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".","noet_",oet_samples,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
  } else {
    saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
  }
}

# print_NjjN_corr_list_strange <- function ( ensemble_id="cA211.30.32", nucleon_type , insertion_type, diagram_type, Gc, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, T ) {
#   nconf <- nconf_f - nconf_i + 1
#   corrlist <- NA
#   for ( i in (nconf_i-1):(nconf_f-1) ) {
#     corr_tmp <- NjjN_corr_strange(
#       file_name = read_file_list( i*4 ),
#       file_dir = read_file_dir( i*4 ),
#       nucleon_type = as.character( nucleon_type ),
#       insertion_type = as.character( insertion_type ),
#       Gc = as.character( Gc ),
#       parity = parity,
#       LT = T )
#     nsrc <- nrow(corr_tmp)
#     corrlist <- c(corrlist, colMeans(corr_tmp))
#   }
#   corrlist <- corrlist[-1]
#   t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
#   
#   if ( Re_Im == "re" ) {
#     corrlist <- cbind(t, remove_im(corrlist))
#   } else if ( Re_Im == "im" ) {
#     corrlist <- cbind(t, remove_re(corrlist))
#   } else if ( Re_Im == "re_im" ) {
#     corr_re <- remove_im(corrlist)
#     corr_im <- remove_re(corrlist)
#     corrlist <- cbind(t, corr_re, corr_im)
#   } else {
#     return("Error: Re_Im must be re, im or re_im")
#   }
#   colnames(corrlist) <- NULL
#   corrlist <- noquote(corrlist)
#   # print(corrlist)
#   saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
# }

# The function print_NjjN_corr_list_niko() was renamed and now is called print_NjjN_corr_list_strange()
print_NjjN_corr_list_strange <- function ( ensemble_id="cA211.30.32", nucleon_type , insertion_type, diagram_type, Gc, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, conf_list = c(), oet_samples = 1, T ) {
  corrlist <- NA
  
  if ( is.null(conf_list) == TRUE ) {
    nconf <- nconf_f - nconf_i + 1
    for ( i in (nconf_i-1):(nconf_f-1) ) {
      corr_tmp <- NjjN_corr_strange(
        file_name = read_file_list( i*4 ),
        file_dir = read_file_dir( i*4 ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  } else {
    nconf <- length(conf_list)
    for ( i in 1:nconf ) {
      corr_tmp <- NjjN_corr_strange(
        file_name = read_file_list( conf_list[i] ),
        file_dir = read_file_dir( conf_list[i] ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  }
  corrlist <- corrlist[-1]
  t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
  
  if ( Re_Im == "re" ) {
    corrlist <- cbind(t, remove_im(corrlist))
  } else if ( Re_Im == "im" ) {
    corrlist <- cbind(t, remove_re(corrlist))
  } else if ( Re_Im == "re_im" ) {
    corr_re <- remove_im(corrlist)
    corr_im <- remove_re(corrlist)
    corrlist <- cbind(t, corr_re, corr_im)
  } else {
    return("Error: Re_Im must be re, im or re_im")
  }
  colnames(corrlist) <- NULL
  corrlist <- noquote(corrlist)
  # print(corrlist)
  if (diagram_type == "W") {
    saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".","noet_",oet_samples,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
  } else {
    saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
  }
}



# corr print function which can process both color crossed (colix = TRUE) and non-color crossed (colix = FALSE) data
print_NjjN_corr_list_light_W <- function ( ensemble_id="cA211.30.32", nucleon_type , insertion_type, diagram_type = "W", Gc, parity = +1, Re_Im = "re", nconf_i = 1, nconf_f = 1, conf_list = c(), oet_samples = 1, colix, single_t, T ) {
  corrlist <- NA
  
  if ( is.null(conf_list) == TRUE ) {
    nconf <- nconf_f - nconf_i + 1
    for ( i in (nconf_i-1):(nconf_f-1) ) {
      corr_tmp <- NjjN_corr_light_W(
        file_name = read_file_list( i*4 ),
        file_dir = read_file_dir( i*4 ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        colix = colix,
        single_t = single_t,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  } else {
    nconf <- length(conf_list)
    for ( i in 1:nconf ) {
      corr_tmp <- NjjN_corr_light_W(
        file_name = read_file_list( conf_list[i] ),
        file_dir = read_file_dir( conf_list[i] ),
        nucleon_type = as.character( nucleon_type ),
        insertion_type = as.character( insertion_type ),
        Gc = as.character( Gc ),
        parity = parity,
        samples = oet_samples,
        colix = colix,
        single_t = single_t,
        LT = T )
      nsrc <- nrow(corr_tmp)
      corrlist <- c(corrlist, colMeans(corr_tmp))
    }
  }
  corrlist <- corrlist[-1]
  t <- rep(c(0:(length(corrlist)/nconf - 1)), nconf)
  
  if ( Re_Im == "re" ) {
    corrlist <- cbind(t, remove_im(corrlist))
  } else if ( Re_Im == "im" ) {
    corrlist <- cbind(t, remove_re(corrlist))
  } else if ( Re_Im == "re_im" ) {
    corr_re <- remove_im(corrlist)
    corr_im <- remove_re(corrlist)
    corrlist <- cbind(t, corr_re, corr_im)
  } else {
    return("Error: Re_Im must be re, im or re_im")
  }
  colnames(corrlist) <- NULL
  corrlist <- noquote(corrlist)

  # Diagram type W:
  saveRDS(corrlist, file = paste0(nucleon_type,"_",insertion_type,"_",nucleon_type,".",ensemble_id,".Nconf_",nconf,".Nsrc_",nsrc,".","noet_",oet_samples,".",diagram_type,".",Gc,".parity_",parity,".",Re_Im,".rds"))
}

average_corr <- function ( corr, ntslices ) {
  options(digits=22)
  corr <- as.numeric(corr[,2])
  ncorr <- length(corr)/ntslices
  corr_sum <- c(rep(0,ntslices))
  corr_av <- c(rep(0,ntslices))
  corr_sd <- c(rep(0,ntslices))
  
  for (i in 1:ntslices) {
    for(j in seq(i,length(corr),ntslices)) {
      corr_sum[i] <- corr_sum[i] + corr[j]
    }
    corr_av[i] <- corr_sum[i]/length(seq(i,length(corr),ntslices))
  }
  
  for (i in 1:ntslices) {
    for(j in seq(i,length(corr),ntslices)) {
      corr_sd[i] <- corr_sd[i] + ( corr[j] - corr_av[i] )^2
    }
    corr_sd[i] <- sqrt((1/length(seq(i,length(corr),ntslices)))*corr_sd[i])
  }
  
  corr_av <- cbind(c(0:(ntslices-1)),corr_sum/ncorr)
  return(cbind(corr_av, corr_sd))
}

fit.const <- function( y,
                       bsamples,
                       x,
                       x.lower, x.upper,
                       start.par,
                       useCov = FALSE,
                       boot.fit = TRUE,
                       fit.method = "optim",
                       autoproceed = FALSE,
                       every,
                       cov_fn = cov,
                       error = sd,
                       ... ) {
  stopifnot(!missing(y))
  stopifnot(!missing(bsamples))
  stopifnot(!missing(x))
  stopifnot(!missing(x.lower))
  stopifnot(!missing(x.upper))
  stopifnot(!missing(start.par))
  
  if (x.upper <= x.lower) {
    stop("Error: x.upper must be > x.lower")
  }
  
  # Prediction function
  prediction_function <- function(par, x, ...) {
    f <- par[1]
    return(f)
  }
  
  # Initial guess for the fit parameter
  initial_guess = function(npar) {
    par <- numeric(npar)
    par <- start.par
    return (par)
  }
  
  fitfun <- prediction_function
  par.guess <- initial_guess(npar = 1)
  
  # Index vector for time slices to be fitted
  i <- c((x.lower+1):(x.upper+1))
  
  args <- list(fn = fitfun,
               par.guess = par.guess,
               y = y[i],
               x = x[i],
               bsamples = bsamples[, i],
               use.minpack.lm = fit.method == 'lm',
               error = error,
               cov_fn = cov_fn,
               x.lower=x.lower, x.upper=x.upper,
               ...)
  
  if (useCov) {
    args$CovMatrix <- cov_fn(bsamples[, i])
  }
  
  res <- do.call(bootstrap.nlsfit, args)
}

rds2cf <- function ( path2rds, Time, symmetrised=FALSE ) {
  cf_tmp <- readRDS( path2rds )
  cf_re <- split(as.double(cf_tmp[,2]), 1:Time)
  cf_im <- split(as.double(cf_tmp[,3]), 1:Time)
  cf_tmp_re <- c()
  cf_tmp_im <- c()
  
  for( i in 1:Time ){
    cf_tmp_re <- cbind(cf_tmp_re, cf_re[[i]])
    cf_tmp_im <- cbind(cf_tmp_im, cf_im[[i]])
  }
  
  cf_tmp <- cf_orig(cf = cf_tmp_re, icf = cf_tmp_im)
  cf_tmp <- cf_meta(nrObs=1, Time=Time, nrStypes = 1, symmetrised=symmetrised, cf_tmp)
  return(cf_tmp)
}

array2cf <- function ( array, Time, parity = +1, symmetrised = FALSE ) {
  cf_tmp_re <- c()
  cf_tmp_im <- c()
  
  if ( parity == +1 ) {
    cf_re <- split( parity*as.double(array[,2]), 1:Time )
    cf_im <- split( parity*as.double(array[,3]), 1:Time )
    tf <- Time
  } else if ( parity == -1 ) {
    cf_re <- split( parity*as.double(rev(array[,2])), 1:Time )
    cf_im <- split( parity*as.double(rev(array[,3])), 1:Time )
    tf <- Time-1
  }
  
  for( i in 1:tf ){
    cf_tmp_re <- cbind(cf_tmp_re, cf_re[[i]])
    cf_tmp_im <- cbind(cf_tmp_im, cf_im[[i]])
  }
  
  if ( parity == -1 ) {
    cf_tmp_re <- cbind(-1*cf_re[[Time]], cf_tmp_re)
    cf_tmp_im <- cbind(-1*cf_im[[Time]], cf_tmp_im)
  }
  
  cf_tmp <- cf_orig(cf = cf_tmp_re, icf = cf_tmp_im)
  cf_tmp <- cf_meta(nrObs=1, Time=Time, nrStypes = 1, symmetrised=symmetrised, cf_tmp)
  return(cf_tmp)
}

convert.corr <- function (corr) {
  corr_tmp <- c()
  for(i in 1:nrow(corr)){
    for(j in 1:ncol(corr)){
      corr_tmp <- rbind(corr_tmp, cbind(j-1, remove_im(corr[i,j]),remove_re(corr[i,j])))
    }
  }
  return(corr_tmp)
}

convert.corr2 <- function (corr) {
  corr_tmp <- c()
  for(i in 1:length(corr[1,,1])){
    for(j in 1:nrow(corr[,i,])){
      for(k in 1:ncol(corr[,i,])){
        corr_tmp <- rbind(corr_tmp, cbind(k-1, remove_im(corr[j,i,k]),remove_re(corr[j,i,k])))
    }
  }
  }
  return(corr_tmp)
}

#plst <- c("-01","00","01")
#momlist <- permutations(n=3,r=3,v=plst,repeats.allowed = TRUE)
#generate_key('p','uu')