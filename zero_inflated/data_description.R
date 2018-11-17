
# This file is used to find some descriptives reported in the paper


#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"



# ------------------------------- #
# data description for zip_data.csv
# ------------------------------- #
zip_data = read.csv(paste0(file_root,"zero_inflated/zip_data.csv"))
# 14 different variables, measuring various defect counts: 
col_zeros = 100*apply(zip_data, 2, function(x) mean(x==0)) 
summary(col_zeros)   # Columnwise zeros account for  ~5% - 99% of values
100*mean(zip_data == 0)     # Total percent is ~60%)

col_nonzero = apply(zip_data, 2, function(x) mean(x[x!=0]))  
summary(col_nonzero) # Mean of non-zero values range from 1.642 - 65.490        (mean range is 1.6 to 65.5)
mean(zip_data[zip_data!=0])        # Global non-zero mean is 13.34443                         (total mean is 13)
median(zip_data[zip_data!=0])      # Global non-zero median is 2                              (total median is 2)



