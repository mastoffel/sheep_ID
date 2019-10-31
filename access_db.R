library(RJDBC)
#### database connection ######
dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

# open connection
con <- DBI::dbConnect(JDBC(driver, driverpath), options)

# names of tables
tbls <- dbGetTables(con)
# names of variables in a table
flds <- dbGetFields(con, "Sheep")
# get a table with all variables
Sheep <- dbGetQuery(con, "Select * from Sheep")