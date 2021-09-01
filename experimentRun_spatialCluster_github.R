library (spatialcluster)

load("simulated_data_list.RData")

simulation_data <- simu_listAndList

xy <- simulation_data[["simu_n"]][["n=100"]][[1]]$xy
data <- simulation_data[["simu_n"]][["n=100"]][[1]]$mat

data_euDistance <- as.matrix(dist(data))
scl <- scl_redcap(xy, data_euDistance, ncl = 2, linkage = "average") #Error!!!!

