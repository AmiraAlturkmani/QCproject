
df_clean <- subset(omim_analysis, select = -c(Symbol, HGNCID, name, avg, std, median))

#omim_analysis is the raw dataset with extra columns than wnated

row.names(df_clean) = omim_analysis$Symbol
