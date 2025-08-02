# TreeLogic.pack

## Introduction

## Installation
The installation process of TreeLogic.pack will depend on "devtools" package.
To avoid some possible issues, it is recommended to start a new Rstudio session and restart R before installing this package.
First we install the "devtools" package using the example code below.
```
install.packages("devtools")
```
Then run:
```
library(devtools)
devtools::install_github("yiwen99/TreeLogic.pack")
library(TreeLogic.pack)
```
And the functions Bottom_up_selection() and Top_down_selection() are ready to be called!

If you want to build the vignettes while installing the package, to ensure successful installation, please do:
(althought the required packages are already included in the Imports section in Description)

```
install.packages("bench")
devtools::install_github("yiwen99/TreeLogic.pack", build_vignettes=TRUE)
library(TreeLogic.pack)
browseVignettes("TreeLogic.pack")
```

If there are messages such as "recommend to update packages xxx", please choose update all. (The installation process may differ for different laptops or different version of R.

## Usage

Please see the details on the usage of the functions in the help pages in this repository. 

To get further assistance on the usage of this package

```
library(TreeLogic.pack)
?Bottom_up_selection
?Top_down_selection
browseVignettes("TreeLogic.pack")
```
