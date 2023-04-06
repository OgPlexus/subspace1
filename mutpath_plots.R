library(tidyverse)
library(igraph, pos='package:base', quietly = T)
library(ggraph, pos='package:base', quietly = T)
library(ungeviz) # From C. Wilke https://github.com/wilkelab/ungeviz
library(corrr)
#### Cleanup functions ----

paths <- function(x){
  if(is.character(x)){
    out <-  switch(x, 
                   "000" = c("001", "010", "100"),
                   "001" = c("011", "101"),
                   "010" = c("011", "110"),
                   "100" = c("101", "110"),
                   "011" = c("111"),
                   "101" = c("111"),
                   "110" = c("111"),
                   character(0)  )
  }else{
    out<- switch(as.character(x), 
                 "0" = c(1, 10, 100),
                 "1" = c(11, 101),
                 "10" = c(11, 110),
                 "100" = c(101, 110),
                 "11" = c(111),
                 "101" = c(111),
                 "110" = c(111),
                 numeric(0)  )
  }
  return(out)
}
bin_allele_count <- function(bin){
  nchar(bin) - nchar(gsub('1', '', bin))
}
bin_int_to_char <- function(x){
  out<- switch(as.character(x), 
               "0" = "000",
               "10" = "010",
               "1" = "001",
               "11" = "011",
               as.character(x) )
  return(out)
}
bin_to_label <- function(x, spp){
  if(spp == "E. coli"){
    out <-  switch(x, 
                   "000" = "PAL",
                   "001" = "PAR",
                   "010" = "PTL",
                   "100" = "LAL",
                   "011" = "PTR",
                   "101" = "LAR",
                   "110" = "LTL",
                   "111" = "LTR",
                   character(0)  )
  } else if(spp =="L. grayi"){
    out <-  switch(x, 
                   "000" = "PAL",
                   "001" = "PAR",
                   "010" = "PTL",
                   "100" = "LAL",
                   "011" = "PTR",
                   "101" = "LAR",
                   "110" = "LTL",
                   "111" = "LTR",
                   character(0)  )
  } else if(spp =="C. muridarum"){
    out <-  switch(x, 
                   "000" = "PEL",
                   "001" = "PER",
                   "010" = "PTL",
                   "100" = "LEL",
                   "011" = "PTR",
                   "101" = "LER",
                   "110" = "LTL",
                   "111" = "LTR",
                   character(0)  )
  } else {out <- "ERROR"}
  out
}

prettylabels <- function(df){
  df |> 
    mutate(species = str_remove(species, " "))|>
    mutate(variable = recode(variable,
                           "IC50" = "IC[50]",
                           "kcat" = "italic('k')[cat]",
                           "ki" ="italic('K')[i]",
                           "km" ="italic('K')[M]",
                           "tm" ="italic(T)[m]"),
           species = recode(species,
                          "L.grayi" = "italic('L. grayi')",
                          "E.coli" = "italic('E. coli')",
                          "C.muridarum" ="italic('C. muridarum')"))
  
}
#### Read data ----

raw_meas <- read.csv("data/measurements.csv")|>
  mutate(haplotype = map_chr(haplotype, bin_int_to_char))|>
  select(-contains("_"), -bisANS)

raw_coef <- read.csv("data/epistasis_coefficients.csv")|>
  mutate(haplotype = map_chr(haplotype, bin_int_to_char))|>
  select(-contains("_"), -bisANS)

raw_epi <- read_csv("data/epistasis_by_order.csv") |>
  select(-contains("_"), -bisANS)|>
  mutate(order = str_remove(pattern = "Order ", order))

#### Mutational path plots ----
mutant_plot <- function(df,xvar='mutcount', yvar='relY', rel=T){
  
  nodelist <- tibble(binary = sort(map_chr(c(0, 1, 10, 100, 11, 101, 110, 111), bin_int_to_char))) %>%
    left_join(., df, by='binary') |>
    mutate(letts = map2_chr(binary, spp, bin_to_label))
  
  xpos <- with(nodelist, get(xvar))
  ypos <- with(nodelist, get(yvar))
  
  templateNet <- tibble(from = nodelist$binary)|>
    mutate(to = map(from, paths) )|>
    unnest(cols=c(to)) 
  
  graf <- igraph::graph_from_data_frame(templateNet, vertices = nodelist$binary)
  
  out <- ggraph::ggraph(graf, layout = data.frame(x = xpos, y = ypos)) +
    ggraph::geom_edge_link0(edge_width=0.45, color='darkcyan', alpha=0.9)+ 
    geom_label(data=nodelist, aes(x=mutcount, y=ypos, label = letts), label.padding = unit(0.09, "lines"), size=2)+
    coord_cartesian(clip='off')+
    theme_minimal()+
    theme(
      panel.grid.minor = element_blank(), 
      panel.grid.major.x  = element_blank(),
      panel.grid.major.y  = element_line(linewidth=0.4, color="grey90"),
      axis.title.y = element_text(angle=0, vjust=0.5),
      axis.ticks.x = element_line(linewidth=0.4, color="grey90"),
      axis.ticks.length.x = unit(-.3, "lines"),
      plot.title = element_text(hjust=0.5)
    )+
    labs(x="Number of mutations", y="")
  if(rel){
    out + 
      scale_y_continuous(breaks = c(0, 1, 2, 3), limits = c(0,3.2)  )+ 
      expand_limits(x=c(-0.1, 3.1))
  } else{
    out
  }
  
}


not_first_col <- function(){
  list(labs(y=""),
       theme(axis.text.y = element_blank())
       )
}
not_first_col_nums <- function(){
  list(labs(y="")
  )
}
not_first_row <- function(){
  list(labs(title="")
  )
}
not_last_row <- function(){
  list(labs(x=""),
    theme(axis.text.x = element_blank(),
          axis.ticks.x =element_blank())
  )
}

processed <- raw_meas |>
  mutate(mutstring = NULL)|>
  pivot_longer(cols = -one_of("species", "haplotype"), values_to = "measure", names_to = "variable")|>
  rename(binary = "haplotype")|>
  mutate(mutcount = map_int(binary,  bin_allele_count),
         spp = species) |>
  group_by(variable, species) |> 
  arrange(binary)|>
  mutate(rel_meas = (measure - min(measure)) / sd(measure)) |>
  ungroup() |>
  arrange(species, variable) |>
  prettylabels()


origami <-   processed|>
  group_by(variable, species) |>
  nest()|>
  mutate(mplot = map(data, mutant_plot, yvar = "rel_meas"))|>
  select(mplot, species, variable)|>
  mutate(mplot = map2(mplot, species, function(m,s) { m + labs(title = parse(text=s)) }),
         mplot = map2(mplot, variable, function(m,v){ m + labs(y = parse(text=v)) }))


origami_panel <- cowplot::plot_grid(
                   origami$mplot[[1]]+not_last_row(), 
                   origami$mplot[[6]]+not_last_row()+not_first_col(), 
                   origami$mplot[[11]]+not_last_row()+not_first_col(), 
                   origami$mplot[[2]]+not_last_row()+not_first_row(), 
                   origami$mplot[[7]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[12]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[3]]+not_last_row()+not_first_row(), 
                   origami$mplot[[8]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[13]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[4]]+not_last_row()+not_first_row(), 
                   origami$mplot[[9]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[14]]+not_last_row()+not_first_col()+not_first_row(), 
                   origami$mplot[[5]]+not_first_row(), 
                   origami$mplot[[10]]+not_first_col()+not_first_row(), 
                   origami$mplot[[15]]+not_first_col()+not_first_row(), 
                   ncol=3, align="v")

ggsave("fig1_mutplot_normalized_meas.pdf", origami_panel, height=10, width = 7.5)  
  


oribase <-   processed|>
  group_by(variable, species) |>
  nest()|>
  mutate(mplot = map(data, mutant_plot, yvar = "measure", rel=F))|>
  select(mplot, species, variable)|>
  mutate(mplot = map2(mplot, species, function(m,s) { m + labs(title = parse(text=s)) }),
         mplot = map2(mplot, variable, function(m,v){ m + labs(y = parse(text=v)) }))


oribase_panel <- cowplot::plot_grid(
  oribase$mplot[[1]]+not_last_row(), 
  oribase$mplot[[6]]+not_last_row()+not_first_col_nums(), 
  oribase$mplot[[11]]+not_last_row()+not_first_col_nums(), 
  oribase$mplot[[2]]+not_last_row()+not_first_row(), 
  oribase$mplot[[7]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[12]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[3]]+not_last_row()+not_first_row(), 
  oribase$mplot[[8]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[13]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[4]]+not_last_row()+not_first_row(), 
  oribase$mplot[[9]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[14]]+not_last_row()+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[5]]+not_first_row(), 
  oribase$mplot[[10]]+not_first_col_nums()+not_first_row(), 
  oribase$mplot[[15]]+not_first_col_nums()+not_first_row(), 
  ncol=3, align="v")

ggsave("altfig1_mutplot_measures.pdf", oribase_panel, height=10, width = 7.5)  

#### Tile plot for measurements ----
bin_to_onelabel <- function(x, spp){
    out <-  switch(x, 
                   "000" = "PAL/PEL",
                   "001" = "PAR/PER",
                   "010" = "PTL",
                   "100" = "LAL/LEL",
                   "011" = "PTR",
                   "101" = "LAR/LER",
                   "110" = "LTL",
                   "111" = "LTR",
                   character(0)  )
  out
}

meas_tileplot <- processed |>
  group_by(variable, species) |> 
  mutate(ranko = min_rank(measure)) |>
  ungroup() |>
  mutate(letts = map2_chr(binary, spp, bin_to_onelabel))|>
  ggplot() + 
  geom_tile(aes(y=species, x=letts, fill = ranko))+
  geom_text(aes(y=species, x=letts, label = ranko), color='white', size=2)+
  facet_wrap(~variable, nrow=1, labeller = label_parsed, scales= 'free_x')+
  scale_fill_viridis(trans = 'reverse')+
  scale_y_discrete(labels = ggplot2:::parse_safe)+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(),
        plot.title.position="plot")+
  labs(y="", x="", fill = "Rank", title = "")

ggsave("figS1_rank_measurement_tiles.pdf", meas_tileplot, height=2, width = 9)

#### Kendall correlation plot ----
cortestp <- function(vec_a, vec_b){
  cor.test(vec_a, vec_b, method='kendall')$p.value
}

facspp <- function(x){factor(x, ordered = T, levels = c("C.muridarum", "E.coli", "L.grayi"), labels = c("italic('C. muridarum')", "italic('E. coli')", "italic('L. grayi')")) }

facmeas <- function(m){ factor(m, ordered = T, levels=c("IC50" ,"kcat" ,"ki", "km", "tm"), labels=c("IC[50]","italic('k')[cat]","italic('K')[i]","italic('K')[M]","italic(T)[m]"))}

prettypairs <- function(df){
  df |>  
    separate(term, into = c("meas_1","species_1"), sep = "_") |>
    separate(term_2, into = c("meas_2","species_2"), sep = "_") |>
    mutate(meas_1 = facmeas(meas_1), meas_2 = facmeas(meas_2))|>
    mutate(species_1 = str_remove(species_1, " "), species_2 = str_remove(species_2, " ")) |>
    mutate(species_1 = facspp(species_1), species_2 = facspp(species_2))
  
}
pairdf <- raw_meas |>
  select(-one_of("mutstring"))|>
  group_by(species) |>
  mutate_at(vars(-group_cols(), "haplotype"),.funs = min_rank) |>
  pivot_wider(id_cols = 'haplotype', names_from = "species", values_from = c(kcat, km, ki, tm, IC50) )|>
  select(-'haplotype')


paircor_k <- pairdf|> 
  colpair_map(cor, method='kendall')|>
  pivot_longer(cols =-c("term"), names_to = "term_2", values_to = "kendall")|>
  prettypairs()

paircor_p <- pairdf|> 
  colpair_map(cortestp)|>
  pivot_longer(cols =-c("term"), names_to = "term_2", values_to = "pval")|>
  prettypairs()|>
  mutate(sig = if_else(pval <0.001, "***", if_else(pval<0.01, "**", if_else(pval<0.05,"*", ""),""),"") ) 

correlplot <- ggplot(paircor_k, aes(x = meas_1, y=meas_2))+
  geom_tile(aes(fill =kendall))+
  geom_text(data = filter(paircor_p, pval <0.05), aes(label=sig))+
  scale_y_discrete(limits=rev, labels = ggplot2:::parse_safe)+
  scale_x_discrete(labels = ggplot2:::parse_safe)+
  scale_fill_gradient2()+
  facet_grid(species_1 ~ species_2, labeller = label_parsed)+
  theme_minimal()+
  theme(strip.text  = element_text(size=10),
        strip.text.y  = element_text(angle=0),
        legend.position = "bottom")+
  labs(x="", y="", fill="Kendall\ncorrelation")

ggsave("altfig3_paircorrelation_rank_plot.pdf", correlplot, height=8, width = 8)

# attempt at reducing
paircor_k_sp <- paircor_k |> filter(as.numeric(species_1) == as.numeric(species_2), meas_1 < meas_2)
paircor_p_sp <- paircor_p |> filter(as.numeric(species_1) == as.numeric(species_2), meas_1 < meas_2)

minforscale <- min(paircor_k$kendall)
maxforscale <- max(paircor_k$kendall)

corr_sppdiag <- ggplot(paircor_k_sp, aes(x = meas_1, y=meas_2))+
  geom_tile(aes(fill =kendall))+
  geom_text(data = filter(paircor_p_sp, pval <0.05), aes(label=sig), size=5)+
  scale_y_discrete(limits=rev, labels = ggplot2:::parse_safe)+
  scale_x_discrete(labels = ggplot2:::parse_safe)+
  scale_fill_gradient2( limits =c(minforscale, maxforscale) )+
  facet_wrap( ~ species_2, labeller=label_parsed)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey90", color=NA),
        panel.grid=element_blank(),
        strip.text  = element_text(size=10, face = "italic"),
        strip.text.y  = element_text(angle=0),
        legend.justification = "right", legend.position = 'bottom',
        legend.title = element_text(size=10, vjust = 0.75))+
  labs(x="", y="", fill="Kendall correlation")



paircor_k_meas <- paircor_k |> filter(as.numeric(species_1) > as.numeric(species_2) , meas_1 == meas_2)
paircor_p_meas <- paircor_p |> filter(as.numeric(species_1) > as.numeric(species_2) , meas_1 == meas_2)

corr_measdiag <- ggplot(paircor_k_meas, aes(x = species_1, y=species_2))+
  geom_tile(aes(fill =kendall))+
  geom_text(data = filter(paircor_p_meas, pval <0.05), aes(label=sig))+
  scale_y_discrete(labels = ggplot2:::parse_safe)+
  scale_x_discrete(limits=rev, labels = ggplot2:::parse_safe)+
  scale_fill_gradient2( limits =c(minforscale, maxforscale) )+
  facet_wrap( ~ meas_2, labeller = label_parsed, nrow=1)+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey90", color=NA),
        panel.grid = element_blank(),
        strip.text  = element_text(),
        strip.text.y  = element_text(angle=0),
        legend.position="none")+
  labs(x="", y="")

alt_kendallplot <- cowplot::plot_grid(corr_sppdiag,
                   corr_measdiag, 
                   ncol=1, rel_heights = c(0.6, 0.3)
                   )

ggsave("fig3_paircorrelation_rank_plot.pdf", alt_kendallplot, height=6.2, width = 8)

 #### Epistasis coefficient plots ----

coefdf <- raw_coef |> 
  pivot_longer(cols = -one_of("species", "haplotype"), values_to = "measure", names_to = "variable")|>
  group_by(variable, species) |> 
  mutate(ranko = min_rank(measure)) |>
  ungroup() |>
  arrange(species, variable) |>
  mutate(haplotype=factor(haplotype, levels = c( "000","001", "010", "100", "011", "101", "110","111"),
                                     labels = c("000", "**1","*1*", "1**","*11","1*1","11*","111"))  )|>
  prettylabels()


rank_line_plot <- ggplot(coefdf, aes(x=haplotype, y=measure, color = species, group = species)) +
  geom_line(linewidth = 0.5)+
  geom_point()+
  facet_wrap(~variable, nrow=1, labeller = label_parsed)+
  scale_x_discrete()+
  scale_color_brewer(labels = ggplot2:::parse_safe, type = "qual")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(linewidth=0.2),
        panel.grid.minor.y = element_line(linewidth=0.2),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(0.5, "lines"),
        plot.title.position="plot",
        legend.position = "top")+
  labs(y="", x="", title = "", color ="")

ggsave("fig4a_epicoef_lines.pdf", rank_line_plot, height=3, width = 10)

rank_tile_plot <- ggplot(coefdf) + 
  geom_tile(aes(y=species, x=haplotype, fill = ranko))+
  geom_text(aes(y=species, x=haplotype, label = ranko), color='white', size=2)+
  facet_wrap(~variable, nrow=1, labeller = label_parsed)+
  scale_fill_viridis(trans = 'reverse', breaks=c(1,3,5,7))+
  scale_y_discrete(labels = ggplot2:::parse_safe)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(),
        plot.title.position="plot",
        legend.position = "none")+
  labs(y="", x="", fill = "Rank", title = "")
ggsave("fig4b_rank_epicoef_tiles.pdf", rank_tile_plot, height=2, width = 10)

#### Epistasis order distribution plots ----

epidf <- raw_epi |>
  pivot_longer(cols = -one_of("species", "order"), values_to = "measure", names_to = "variable")|>
  prettylabels()

orderplot <- ggplot(epidf)+geom_col(aes(x=order, y=measure, fill=order))+
  facet_grid(variable~species, switch='y', scales="free_y", labeller = label_parsed)+
  theme_minimal()+
  scale_fill_ordinal()+
  scale_y_continuous(breaks=c(0,0.2, 0.4, 0.6, 0.8,1))+
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_line(linewidth=0.4),
        strip.placement = "outside",
        strip.text.y.left  = element_text(angle=0),
        panel.spacing = unit(2, "lines"),
        legend.position="none",
        plot.title.position = "plot")+
  labs(y="", x="Order of coefficient", title="")


ggsave("fig5_epi_by_order.pdf", orderplot, height=7, width = 4.5)
