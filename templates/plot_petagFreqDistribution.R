#! /bin/env Rscript --no-save

library(tidyverse)

file.path('$tag_directory_path', 'petag.FreqDistribution_1000.txt') %>%
	read.delim(header=TRUE) %>%
	set_names(c('DISTANCE','FREQUENCY')) %>%
	head(n=-1) %>% # remove the 'More than...' row
	subset(FREQUENCY>0) %>%
	mutate_all(as.numeric) %>%
	ggplot() +
	aes(x=DISTANCE, y=FREQUENCY) +
	geom_line(alpha=0.5) +
	geom_smooth(fill=NA) +
	scale_x_continuous(trans='log10') +
	scale_y_continuous(trans='log10', labels=scales::scientific) +
	labs(x='Genomic distance (bp)', y='Interaction frequency') +
	theme_bw() +
	theme(aspect.ratio=1) -> gg

ggsave(filename=sprintf('%s.%s', '$dataset_name', '$output_format'), plot=gg)
