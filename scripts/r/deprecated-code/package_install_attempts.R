# Package installation for working with LandsatTS
# Calum Hoad

# Things I attempted to get install working -------
# install.packages("devtools")
devtools::install_github("logan-berner/LandsatTS", build_vignettes = TRUE)

#devtools::install_github('r-lib/systemfonts')
#install.packages("systemfonts", dependencies = TRUE)
#library(systemfonts)
#tinytex::reinstall_tinytex()  # can't find tlmgr

#devtools::install_local("C:/Users/s1437405/Documents/PhD_Local/c1-analyses/TinyTeX-1-v2023.10.zip")
#tinytex::install_tinytex()

# install.packages('tinytex') < ----- this one worked
library(tinytex)
install.packages('pdflatex')

Sys.which('pdflatex')
.libPaths()

Sys.getenv('PATH')

Sys.setenv(PATH=paste(Sys.getenv("PATH"),"C:/Users/s1437405/AppData/Local/Programs/MiKTeX2/miktex/bin/x64",sep=";"))

Sys.which('pdflatex')

Sys.getenv('PATH')

Sys.getenv('PATH')

library(LandsatTS)

#ee_install_set_pyenv('C:/Users/s1437405/Anaconda3')