# Circos-Plot-of-cancer-data
# Libraries used to create a circos plot 
library("RCircos")
library(circlize)
library(dplyr)
library(biomaRt) 

# Input is an ideal UCSC.HG38.Human.CytoBandIdeogram
data("UCSC.HG38.Human.CytoBandIdeogram")

cyto.info = UCSC.HG38.Human.CytoBandIdeogram

RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()
![ideogram](https://user-images.githubusercontent.com/110582335/197979302-187ee3e5-00e5-4b01-961f-c698e99e40dd.png)

# Modifying the outer band
cyto.info = UCSC.HG38.Human.CytoBandIdeogram

cyto.info$Name = NA

cyto.info$Stain = NA

RCircos.Set.Core.Components(cyto.info, 
                            chr.exclude=NULL, 
                            tracks.inside=10, 
                            tracks.outside=0)

# Modifying a circos plot according to user demand 
ideo = RCircos.Get.Plot.Ideogram()

ideo$BandColor = 'salmon'

num = which(ideo$Chromosome == 'chrX')

ideo[num, 'BandColor'] = 'Orange'

num = which(ideo$Chromosome == 'chrY')

ideo[num, 'BandColor'] = 'Yellow'


RCircos.Reset.Plot.Ideogram(ideo)

RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()
![id](https://user-images.githubusercontent.com/110582335/197979781-2bbc2ef2-d295-43dd-befb-6a94a8d57381.png)


# Input of dataframe created in Differential Gene expression 
mat = readRDS("mat_g.rds")

m = useMart('ensembl', dataset='hsapiens_gene_ensembl')

# Create co-ordinates 
coords = getBM(attributes=c('chromosome_name', 'start_position', 
                            'end_position', 'hgnc_symbol'),
               filters = c('hgnc_symbol'),
               values = list(rownames(mat)),
               mart = m)

View(coords)
![image](https://user-images.githubusercontent.com/110582335/200104864-c3e048d4-a40a-46e7-9068-8e4cdcaf1922.png)

write.csv(coords, file = 'coords.csv')

# Assigning Co-ordinates to chromosome names
coords$chromosome_name = paste0('chr', coords$chromosome_name)

chr_order = unique(cyto.info$Chromosome)

coords$chromosome_name = factor(coords$chromosome_name, levels = chr_order)

# Further modification
num = which(is.na(coords$chromosome_name))

coords = coords[-num, ]

up = which((mat$pval < 0.01) &
             (mat$log2FC > 1))
upmat = mat[up, ]

num = which(coords$hgnc_symbol %in% rownames(upmat))

coords1 = coords[num, ]

table(coords1[,1])
![image](https://user-images.githubusercontent.com/110582335/200104844-99b8ffa1-6822-4881-9139-d1252eea7047.png)

RCircos.Gene.Name.Plot(coords1, name.col=4, track.num = 2, side = "in",
                       is.sorted = F)
     ![image](https://user-images.githubusercontent.com/110582335/200104751-f3f34957-d06a-495e-b5a0-ffabf6605057.png)
                  

# to view the number of labels for each chromosome.
genes = intersect(rownames(mat), coords$hgnc_symbol)

mat1 = mat[genes, ]

# Create a data frame to use to plot a final circos-plot 
df = cbind.data.frame(rownames(mat1), mat1[, c(1,2,4)])
colnames(df)[1] = 'hgnc_symbol'

data = merge(coords, df, by = 'hgnc_symbol', all.x = T)

data = data[, c('chromosome_name', 'start_position',
                'end_position', 'hgnc_symbol',
                'meanTumor', 'meanControl', 'log2FC')]
                
# Plot a circos-plot 
RCircos.Heatmap.Plot(data, data.col = 7, track.num = 6, side = "in",
                     min.value = -3, max.value = 6, 
                     is.sorted = F)

RC.param = RCircos.Get.Plot.Parameters()

RCircos.Reset.Plot.Parameters(RC.param)
![image](https://user-images.githubusercontent.com/110582335/200104612-5d46e401-f2fd-416c-8b11-245c531c92c8.png)


# Interpretation :
Outer circle indicates number of chromosome pairs 
2nd circle indicates presence of genes on respective chromosomes 
Inner circle represent a heatmap whereas blank spaces in it shows absence of genes .
