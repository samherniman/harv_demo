#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| eval: false
install.packages(c("pak", "neonOS", "neonUtilities", "here"))
pak::pak("samherniman/easytrees")
#
#
#
#
#
#| warning: false
#| error: false
future::plan(future::multisession, workers = 8L)
lidR::set_lidr_threads(4)
lasR::set_parallel_strategy(lasR::nested(2, 4L))

output_dir <- here::here("data/derivative")
#
#
#
#
#
#| eval: false

neonUtilities::byTileAOP(
  dpID = "DP1.30003.001",
  site = "HARV",
  year = "2022",
  easting = 732183,
  northing = 4713265,
  buffer = 200,
  savepath = here::here("data/raw/las/harv")
)
#
#
#
#
#
ctg <- lidR::readLAScatalog(here::here("data/raw/las/harv"), recursive = TRUE)
# plot(ctg, mapview = TRUE)
ctg
#
#
#
#
#
#| output: false
#| warning: false
#| error: false
pipeline <- 
  easytrees::setup_pipeline(output_dir) |> 
  easytrees::create_pipeline()

ans = lasR::exec(pipeline, on = ctg)
#
#
#
#
#
#
#
#| output: false
#| warning: false
#| error: false
ctg <- lidR::readLAScatalog(
  fs::path(output_dir, "normalized")
)

chm <- terra::rast(fs::path(output_dir, "chm", "chm.tif"))

f <- function(x) {
  y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3
  y[x < 2] <- 3
  y[x > 20] <- 5
  return(y)
}
lidR::opt_chunk_size(ctg) <- 5000
#
#
#
ctg
#
#
#
#
#
#
#| output: false
#| warning: false
#| error: false
ttops <- lidR::locate_trees(ctg, lidR::lmf(f), uniqueness = "bitmerge")
ttops$treeID <- 1:nrow(ttops)
sf::st_write(ttops, fs::path(output_dir, "treetops", "treetops.fgb"), append = FALSE)
#
#
#
ttops
#
#
#
#
#
#
#| output: false
#| warning: false
#| error: false
#| eval: false
algo3 <- lidR::watershed(chm, th_tree = 2, tol = 1, ext = 1)
lidR::opt_laz_compression(ctg) <- TRUE
lidR::opt_output_files(ctg) <- paste0(output_dir, "/clouds/water_{XCENTER}_{YCENTER}_d")
ctg <- lidR::segment_trees(ctg, algo3, attribute = "IDwater")
#
#
#
ctg
#
#
#
#
#
#
#
ctg <- lidR::readLAScatalog(
  fs::path(output_dir, "clouds")
)
lidR::opt_chunk_size(ctg) <- 50000

#
#
#
library(lidR)
# las <- lidR::readLAS(ctg[1,])
cm_las <- lidR::crown_metrics(ctg, func = .stdtreemetrics, attribute = "IDwater", geom = "concave") |> 
  dplyr::filter(convhull_area < 400) |> 
  sf::st_make_valid()
# removed all trees that are greater than 400 square meters in area
# because those ones usually have really weird hulls
#
#
#
sf::st_write(cm_las, fs::path(output_dir, "hulls.fgb"), append = FALSE)
#
#
#
#
mapview::mapview(cm_las)
plot(
  cm_las,
  #  "IDwater", 
   col = pastel.colors(200),
   xlim = c(731001, 731500),
   ylim = c(4713000, 4713500)
   )
#
#
#
