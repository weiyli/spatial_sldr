ed_table1 <- data.frame(
  Event = c(
    "Hurricane Harvey",
    "Hurricane Dorian",
    "Winter Storm Uri",
    "Kincade Fire",
    "COVID-19"
  ),
  `Baseline window` = c(
    "18-24 August 2017",
    "26 August-1 September 2019",
    "30 January-9 February 2021",
    "8-22 October 2019",
    "10-16 February 2020"
  ),
  `Event window` = c(
    "25 August-1 September 2017",
    "2-5 September 2019",
    "10-20 February 2021",
    "23 October-6 November 2019",
    "30 March-5 April 2020"
  ),
  `Study region` = c(
    "Houston, Texas",
    "Jacksonville, Florida",
    "Houston, Texas",
    "Santa Rosa-Petaluma, California",
    "12 U.S. metropolitan areas"
  ),
  check.names = FALSE
)

write.csv(
  ed_table1,
  "D:/ood/Data/spatial_sldr/sldr_ED_Table1.csv",
  row.names = FALSE
)


