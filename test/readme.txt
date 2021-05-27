Descriptions of test files:

test-biglasso-no-cv.R
    Compare the penalized Cox regressions implementated in biglasso::biglasso() 
    and glmnet::glmnet() in order to ensure that I know how to use biglasso.
    Note that this test is not leveraging the utility offered by 
    bigmemory::big.matrix since the test data is fully loaded into R (this is
    fine since it wasn't the point of this test).
    
test-biglasso-cv.R
    Compare the output of our new cross-validation function (cv.biglasso.cox()) 
    with the values obtained by glmnet (cv.glmnet()) in a memory-efficient way
    (i.e. taking full advtange of bigmemory::big.matrix objects).

test-biglasso-cv-old.R
    Compare the output of our new cross-validation function (cv.biglasso.cox()) 
    with the values obtained by glmnet (cv.glmnet()).
    OLD VERSION: Does not leverage the utility offered by bigmemory::big.matrix.

structure-data-for-biglasso.R
    Build a framework to take the data in its original format and convert
    it into a format more directly ammenable to biglasso::biglasso(),
        1) (optional/temporary) deal with NA/missing values,
        2) convert the raw predictor columns into the corresponding design matrix,
    where both steps ought to be done in a way that is memory efficient (since
    the full data is at least 118000 rows * 58000 columns * 8 bytes per entry
    = 55 GB)

test-big-matrix.R
    Understand how to use bigmemory::big.matrix in order to have a memory
    efficient data structure to use with biglasso::biglasso() (i.e.
    reading the data via read.big.matrix to avoid having the data enter
    memory as it would if it was read with read.csv, data.table::fread, etc.)

test-big-matrix-cv-subsetting.R 
    Figuring out how to generate training/test subsets of the data while 
    keeping the data (and its subsets) big.matrix objects for use in a
    biglasso::biglasso() cross-validation function (i.e. without 
    converting the data/subsets to a data.frame or matrix as an intermediate 
    step prior to re-converting back to a big.matrix).
    
test-old.R
    An older test file whose framework may have something useful for me in the
    future.
    
    
    
    
    
    