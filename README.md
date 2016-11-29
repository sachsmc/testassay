[![Travis-CI Build Status](https://travis-ci.org/sachsmc/testassay.svg?branch=master)](https://travis-ci.org/sachsmc/testassay)

# testassay

An R package to facilitate assay validation in a hypothesis testing framework. This is a companion to the following paper, currently under review:

Michael P Fay, Michael C Sachs, and Kazutoyo Miura. *A Hypothesis Testing Framework for Validating an Assay for Precision*. (2016), Under Review.

### Abstract

A common way of validating a biological assay for is through a procedure, where m levels of an analyte are measured with n replicates at each level, and if all m estimates of the coefficient of variation (CV) are less than some prespecified level, then the assay is declared validated for precision within the range of the m analyte levels. Two limitations of this procedure are: there is no clear statistical statement of precision upon passing, and it is unclear how to modify the procedure for assays with constant standard deviation. We provide tools to convert such a procedure into a set of m hypothesis tests. This reframing motivates the m:n:q procedure, which upon completion delivers a 100q% upper confidence limit on the CV. Additionally, for a post-validation assay output of y, the method gives an ``effective standard deviation interval'' of log(y) plus or minus r, which is a 68% confidence interval on log(mu), where mu is the expected value of the assay output for that sample. Further, the m:n:q procedure can be straightfowardly applied to constant standard deviation assays. We illustrate these tools by applying them to a growth inhibition assay.

### Details

Read the vignette to see how the package is used <https://sachsmc.github.io/testassay>. You are welcome to send us feedback using [Github issues](https://github.com/sachsmc/testassay/issues).
