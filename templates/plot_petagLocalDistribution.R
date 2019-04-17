#! /bin/env Rscript --no-save

library(tidyverse)

file.path('$tag_directory_path', 'petag.LocalDistribution.txt') %>%
	read.delim(header=TRUE) %>%
	set_names(c('DISTANCE','Same','Opposite')) %>%
	gather(key='STRAND', value='FREQUENCY', -DISTANCE) %>%
	subset(FREQUENCY>0) %>%
	ggplot() +
	aes(x=DISTANCE, y=FREQUENCY, colour=STRAND) +
	geom_line(alpha=0.5) +
	geom_smooth(fill=NA) +
	labs(x='Genomic distance (bp)', y='Frequency (m)', colour='Stands') +
	coord_cartesian(xlim=5e3*c(-1,1)) +
	theme_bw() +
	theme(aspect.ratio=1) -> gg

ggsave(filename=sprintf('%s.%s', '$dataset_name', '$output_format'), plot=gg)
