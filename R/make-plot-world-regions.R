
# NB: I use UNAIDS classifications for MENA
# I use UN bloc defns for Western Europe (which includes Greece)

# plot world regions
infile.georeg <- file.path('/Users/alexb/Box Sync/Roadmap','misc','NEWGEO_220713.csv')
geo <- data.table(read.csv(infile.georeg))
geo[geo$Alpha_2_code %in%c('AM','AZ','BY','GE','MD','RU','UA'),WRLD:='EEurope']
geo[geo$WRLD=='Oceania',WRLD:='Other']
geo[geo$Alpha_2_code %in%c('AU','NZ'),WRLD:='Oceania']
setnames(geo,c('CNTRY','WRLD'),c('CNTRY_born','WRLD_born'))

geo[, GEOREG:="Other"]
geo[WRLD_born %in% c("WEurope","NorthAm","Oceania"), GEOREG:="W.Europe,\nN.America,Oceania"]
geo[WRLD_born %in% c("EEurope", "CEurope"), GEOREG:="E. & C. Europe"]
geo[WRLD_born %in% c("LaAmCar"), GEOREG:="S. America &\nCaribbean"]
geo[WRLD_born %in% c("DutchCarSuriname"), GEOREG:="Suriname &\nDutch Caribbean"]
geo[WRLD_born %in% c("MENA"), GEOREG:="MENA"]
geo[Alpha_2_code=="NL", GEOREG:="Netherlands"]

# Retrievethe map data
#map <- data.table(map_data("world"))
# map data with alpha2 codes from https://github.com/Thom-J-H/map_Gap_2_Tidy/blob/main/world_map2.rds
# https://rpubs.com/Thom_JH/798825
map <- data.table(readRDS('/Users/alexb/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Roadmap/sources/ethnicity_analysis/paper2_figures/world_map2.rds'))

# Compute the centroid as the mean longitude and lattitude
# Used as label coordinate for country's names
region.lab.data <- map[, list(long=mean(long),
                              lat=mean(lat)),
                              by='code_2']
region.lab.data <- merge(region.lab.data,geo,by.x='code_2',by.y='Alpha_2_code',all.x=T)
region.lab.data <- subset(region.lab.data,!is.na(code_2))


map <- merge(map,geo,by.x='code_2',by.y='Alpha_2_code',all.x=T)
map <- subset(map,!is.na(code_2))
map[country=='Canary Islands', GEOREG:= 'W.Europe,\nN.America,Oceania']
map[country=='Kosovo', GEOREG:= 'E. & C. Europe']
map[country=='Namibia', GEOREG:= 'Other']
map[country=='Tanzania', GEOREG:= 'Other']

map[, GEOREG:= factor(GEOREG,levels=c('Netherlands','W.Europe,\nN.America,Oceania','Suriname &\nDutch Caribbean',
                                                                         'S. America &\nCaribbean','E. & C. Europe','MENA','Other'))]

g <- ggplot(subset(map,code_2!='AQ'), aes(x = long, y = lat)) +
  geom_polygon(aes( group = group, fill = GEOREG))+
  #geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  scale_fill_npg() +
  labs(fill='') +
  theme_void()+
  theme(legend.position = "right")
ggsave(file = paste0(outfile.base,'-map_georegs.pdf'), g, w = 8, h = 4)
ggsave(file = paste0(outfile.base,'-map_georegs.png'), g, w = 8, h = 4)

pal <- pal_npg('nrc')(7)[c(3,4)]
#map_crop <- st_crop(map, xmin = -73, xmax = -51,
#                          ymin = -0.29, ymax = 15)
g <- ggplot(subset(map,code_2!='AQ' & GEOREG %in% c('Suriname &\nDutch Caribbean','S. America &\nCaribbean') &
                     long > -73 & long < -51 & lat > -0.29 & lat < 15), aes(x = long, y = lat)) +
  geom_polygon(aes( group = group, fill = GEOREG))+
  #geom_text(aes(label = region), data = region.lab.data,  size = 3, hjust = 0.5)+
  #scale_fill_viridis_d()+
  #scale_fill_npg() +
  scale_fill_manual(values=pal) +
  labs(fill='') +
  theme_void()+
  theme(legend.position = "right")
ggsave(file = paste0(outfile.base,'-map_georegs_Suriname_DC.pdf'), g, w = 8, h = 4)
ggsave(file = paste0(outfile.base,'-map_georegs_Suriname_DC.png'), g, w = 8, h = 4)

