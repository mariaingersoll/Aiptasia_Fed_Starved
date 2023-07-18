#Beginning from line 551 from "Aiptasia_FS_2022_clean.Rmd"
#GO Analysis
#- this analysis is performed in GO_MWU.R script

#-read back in your go files, these will go into the dendrogram analysis
#- make sure your files are in unix first, so open and save them in BBEdit, and rename them with the proper "isogroup" instead of aipgene
#- load all the necessary files into your working directory

setwd("/Users/mariaingersoll/Desktop/BU_Research/Davies-Gilmore/Aip_FS_2022/Aiptasia_Fed_Starved/Data_Files/") 

#APO FIRST
    # CC
input="apo_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#208 terms at 10% FDR
resultsCC=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as CCapo_012223.pdf as portrait with 30x10 dimensions in Figures
resultsCC


    # BP
input="apo_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#796
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as BPapo_012223 as portrait with 100x10 dimensions in Figures
resultsBP

  # MF
input="apo_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#276
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as MFapo_012223.pdf as portrait with 45x10 dimenstions in Figures
resultsMF

#Now Sym
  # CC
input="sym_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="CC"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#151
resultsCC=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as CCsym_012223.pdf as portrait with 30x10 dimensions in Figures
resultsCC

  # BP
input="sym_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="BP"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#488
resultsBP=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as BPsym_012223.pdf as portrait with 76x10 dimensions in Figures
resultsBP

  # MF
input="sym_GO_22.csv"
goAnnotations="aiptasia_iso2go.tab"
goDatabase="go.obo"
goDivision="MF"
source("gomwu.functions.R")

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl",
           largest=0.1,
           smallest=5,
           clusterCutHeight=0.25)
#173
resultsMF=gomwuPlot(input,goAnnotations,goDivision,
                    absValue=1,
                    level1=0.1,
                    level2=0.05,
                    level3=0.01,
                    txtsize=1.5,
                    treeHeight=0.5,
)
#saved as MFsym_012223 as portrait with dimensions 25x10 in Figures
resultsMF

###Now return to line 559 of "Aiptasia_FS_2022_clean.Rmd"