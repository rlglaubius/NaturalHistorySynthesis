library(first90)

## Set of countries to use in simulation
metadata = data.frame(
  country = c('Cameroon', "Cote d'Ivoire", 'Malawi', 'Namibia', 'Nigeria', 'Swaziland', 'Tanzania', 'Uganda', 'Zambia', 'Zimbabwe'),
  pjnz = c('Cameroon2019.PJNZ',
           'CotedIvoire2019.PJNZ',
           'Malawi2019.PJNZ',
           'Nigeria2019.PJNZ',
           'Namibia2019.PJNZ',
           'Eswatini2019.PJNZ',
           'Tanzania2019.PJNZ',
           'Uganda2019.pjnz',
           'Zambia2019.PJNZ',
           'Zimbabwe2019.PJNZ'),
  alpha = c('CMR', 'CIV', 'MWI', 'NAM', 'NGA', 'SWZ', 'TZA', 'UGA', 'ZMB', 'ZWE'),
  stringsAsFactors=FALSE)

pjnz.path = 'Data/PJNZ'

spec.info = lapply(sprintf('%s/%s', pjnz.path, metadata$pjnz), prepare_inputs)

names(spec.info) = metadata$alpha

rm(pjnz.path)
