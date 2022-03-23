.onLoad <- function(lib, pkg) {

}

.onAttach <- function(libname, pkgname) {
    db_example_track_db <- sprintf("%s/%s", system.file("trackdb", package = "shaman"), "test")
    if (!file.exists(db_example_track_db)) {
        message("building example db...")
        current_wd <- getwd()
        setwd(system.file("", package = "shaman"))
        system(paste("tar -zxf trackdb.tar.gz"))
        setwd(current_wd)
    }
    .shaman_load_config(system.file("config", package = "shaman"))
}
