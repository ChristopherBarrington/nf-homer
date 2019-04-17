#! /bin/env Rscript --no-save

library(tidyverse)

file.path('$tag_directory_path', sprintf('petagRestrictionDistribution.%s.mis%s.txt', '$restriction_site', '$mismatches')) %>%
	read.delim(header=TRUE, stringsAsFactors=FALSE) %>%
	set_names(c('DISTANCE','Forward','Reverse')) %>%
	gather(key='STRAND', value='FREQUENCY', -DISTANCE) %>%
	subset(FREQUENCY>0) %>%
	(function(x) {as_tibble(x) %>% print(); x}) %>%
	ggplot() +
	aes(x=DISTANCE, y=FREQUENCY, colour=STRAND) +
	geom_line() +
	labs(x=sprintf('Genomic distance to %s (bp)', '$restriction_site'), y='Read counts', colour='Strand') +
	theme_bw() +
	theme(aspect.ratio=1, legend.title=element_blank()) -> gg

ggsave(filename=sprintf('%s.%s', '$dataset_name', '$output_format'), plot=gg)
