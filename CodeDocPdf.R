pack <- "MsdeParEst"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))

R CMD Rd2pdf "~/Montages_reseau/home_maiage/Recherche/MsdeParEst/MsdeParEst/MsdeParEst"
R CMD Rd2pdf "/Users/mdelattre/Documents/R/MsdeParEst/MsdeParEst/man"
