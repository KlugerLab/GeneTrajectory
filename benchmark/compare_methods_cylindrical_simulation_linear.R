source("./order_genes_functions_linear.R")

##### Prepare the hyperparameter list
param_list <- list()
for (N_cells in c(500, 1000, 2500)){
  for (sparsity in c(0.025, 0.05, 0.1, 0.2)){
    for (i in 1:10){
      param_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]] <- c(N_cells, sparsity, i)
    }
  }
}

##### Load preprocessed data
data_S_list <- readRDS("./data/data_S_list_cylinder.rds")

##### Order genes based on the cell pseudotime inferred by different methods
monocle3_linear_GT_cor <- lapply(param_list, compare_GT_Monocle3_linear)
saveRDS(monocle3_linear_GT_cor, "monocle3_linear_GT_cor.rds")

monocle2_linear_GT_cor <- lapply(param_list, compare_GT_Monocle2_linear)
saveRDS(monocle2_linear_GT_cor, "monocle2_linear_GT_cor.rds")

PAGA_linear_GT_cor <- lapply(param_list, compare_GT_PAGA_linear)
saveRDS(PAGA_linear_GT_cor, "PAGA_linear_GT_cor.rds")

GeneTrajectory_linear_GT_cor <- lapply(param_list, compare_GT_GeneTrajectory_linear)
saveRDS(GeneTrajectory_linear_GT_cor, "GeneTrajectory_linear_GT_cor.rds")

SlingShot_linear_GT_cor <- lapply(param_list, compare_GT_SlingShot_linear)
saveRDS(SlingShot_linear_GT_cor, "SlingShot_linear_GT_cor.rds")

##### Monocle3
monocle3_linear_GT_table <- data.frame(name = names(monocle3_linear_GT_cor),
                                       GT = unlist(monocle3_linear_GT_cor))
monocle3_linear_GT_table$NC <- unlist(strsplit(monocle3_linear_GT_table$name, "_"))[c(1:nrow(monocle3_linear_GT_table))*6-4]
monocle3_linear_GT_table$NC <- factor(as.character(monocle3_linear_GT_table$NC), levels = c('500', '1000', '2500'))
monocle3_linear_GT_table$SPARSITY <- unlist(strsplit(monocle3_linear_GT_table$name, "_"))[c(1:nrow(monocle3_linear_GT_table))*6-2]
monocle3_linear_GT_table$replicate <- unlist(strsplit(monocle3_linear_GT_table$name, "_"))[c(1:nrow(monocle3_linear_GT_table))*6]


##### Monocle2
monocle2_linear_GT_table <- data.frame(name = names(monocle2_linear_GT_cor),
                                       GT = unlist(monocle2_linear_GT_cor))
monocle2_linear_GT_table$NC <- unlist(strsplit(monocle2_linear_GT_table$name, "_"))[c(1:nrow(monocle2_linear_GT_table))*6-4]
monocle2_linear_GT_table$NC <- factor(as.character(monocle2_linear_GT_table$NC), levels = c('500', '1000', '2500'))
monocle2_linear_GT_table$SPARSITY <- unlist(strsplit(monocle2_linear_GT_table$name, "_"))[c(1:nrow(monocle2_linear_GT_table))*6-2]
monocle2_linear_GT_table$replicate <- unlist(strsplit(monocle2_linear_GT_table$name, "_"))[c(1:nrow(monocle2_linear_GT_table))*6]


##### PAGA|DPT
PAGA_linear_GT_table <- data.frame(name = names(PAGA_linear_GT_cor),
                                   GT = unlist(PAGA_linear_GT_cor))
PAGA_linear_GT_table$NC <- unlist(strsplit(PAGA_linear_GT_table$name, "_"))[c(1:nrow(PAGA_linear_GT_table))*6-4]
PAGA_linear_GT_table$NC <- factor(as.character(PAGA_linear_GT_table$NC), levels = c('500', '1000', '2500'))
PAGA_linear_GT_table$SPARSITY <- unlist(strsplit(PAGA_linear_GT_table$name, "_"))[c(1:nrow(PAGA_linear_GT_table))*6-2]
PAGA_linear_GT_table$replicate <- unlist(strsplit(PAGA_linear_GT_table$name, "_"))[c(1:nrow(PAGA_linear_GT_table))*6]



##### GeneTrajectory
GeneTrajectory_linear_GT_table <- data.frame(name = names(GeneTrajectory_linear_GT_cor),
                                             GT = unlist(GeneTrajectory_linear_GT_cor))
GeneTrajectory_linear_GT_table$NC <- unlist(strsplit(GeneTrajectory_linear_GT_table$name, "_"))[c(1:nrow(GeneTrajectory_linear_GT_table))*6-4]
GeneTrajectory_linear_GT_table$NC <- factor(as.character(GeneTrajectory_linear_GT_table$NC), levels = c('500', '1000', '2500'))
GeneTrajectory_linear_GT_table$SPARSITY <- unlist(strsplit(GeneTrajectory_linear_GT_table$name, "_"))[c(1:nrow(GeneTrajectory_linear_GT_table))*6-2]
GeneTrajectory_linear_GT_table$replicate <- unlist(strsplit(GeneTrajectory_linear_GT_table$name, "_"))[c(1:nrow(GeneTrajectory_linear_GT_table))*6]



##### CellRank
CellRank_linear_GT_cor <- readRDS("./data/benchmark/CellRank_linear_GT_cor.rds")
CellRank_linear_GT_table <- data.frame(name = names(CellRank_linear_GT_cor),
                                       GT = unlist(CellRank_linear_GT_cor))
CellRank_linear_GT_table$NC <- unlist(strsplit(CellRank_linear_GT_table$name, "_"))[c(1:nrow(CellRank_linear_GT_table))*6-4]
CellRank_linear_GT_table$NC <- factor(as.character(CellRank_linear_GT_table$NC), levels = c('500', '1000', '2500'))
CellRank_linear_GT_table$SPARSITY <- unlist(strsplit(CellRank_linear_GT_table$name, "_"))[c(1:nrow(CellRank_linear_GT_table))*6-2]
CellRank_linear_GT_table$replicate <- unlist(strsplit(CellRank_linear_GT_table$name, "_"))[c(1:nrow(CellRank_linear_GT_table))*6]


##### SlingShot
SlingShot_linear_GT_cor <- readRDS("/banach1/rq25/GeneTrajectory_data/paper/SlingShot_linear_GT_cor.rds")
SlingShot_linear_GT_table <- data.frame(name = names(SlingShot_linear_GT_cor),
                                        GT = unlist(SlingShot_linear_GT_cor))
SlingShot_linear_GT_table$NC <- unlist(strsplit(SlingShot_linear_GT_table$name, "_"))[c(1:nrow(SlingShot_linear_GT_table))*6-4]
SlingShot_linear_GT_table$NC <- factor(as.character(SlingShot_linear_GT_table$NC), levels = c('500', '1000', '2500'))
SlingShot_linear_GT_table$SPARSITY <- unlist(strsplit(SlingShot_linear_GT_table$name, "_"))[c(1:nrow(SlingShot_linear_GT_table))*6-2]
SlingShot_linear_GT_table$replicate <- unlist(strsplit(SlingShot_linear_GT_table$name, "_"))[c(1:nrow(SlingShot_linear_GT_table))*6]


##### Combine all tables
GeneTrajectory_linear_GT_table$method <- "GeneTrajectory"
PAGA_linear_GT_table$method <- "PAGA|DPT"
CellRank_linear_GT_table$method <- "CellRank"
monocle3_linear_GT_table$method <- "Monocle3"
monocle2_linear_GT_table$method <- "Monocle2"
SlingShot_linear_GT_table$method <- "SlingShot"

all_tables <- list(GeneTrajectory_linear_GT_table,
                   PAGA_linear_GT_table,
                   CellRank_linear_GT_table,
                   monocle3_linear_GT_table,
                   monocle2_linear_GT_table,
                   SlingShot_linear_GT_table)
combined_table <- do.call(rbind, all_tables)
combined_table$method <- factor(combined_table$method, levels = c("GeneTrajectory", "PAGA|DPT", "SlingShot", "Monocle2", "Monocle3", "CellRank"))
ggplot(data = combined_table, aes(x = method, y = GT, fill = method)) + geom_boxplot() + facet_grid(NC ~ SPARSITY, margins = FALSE) & RotatedAxis()


data_new <- combined_table %>%
  group_by(method, SPARSITY, NC) %>%
  summarise( 
    n=n(),
    mean=mean(GT),
    sd=sd(GT)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))


ggplot(data = data_new) + 
  facet_grid(SPARSITY ~ NC, margins = FALSE) +
  geom_bar(aes(x = method, y = mean, fill = method), stat="identity", alpha=1, width = .75) +
  geom_errorbar(aes(x=method, ymin=mean-se, ymax=mean+se), width=0.5, colour="black", alpha=0.9, size=.25) +
  coord_flip() + scale_x_discrete(limits=rev) + ylab("Spearman correlation")




