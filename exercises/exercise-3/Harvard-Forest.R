## ----simulate_data------------------------------------------------------------
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)
n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))

x <- cbind(1, rnorm(n))

B <- as.matrix(c(1,5))

sigma.sq <- 5
tau.sq <- 1
phi <- 3/0.5

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, x%*%B + w, sqrt(tau.sq))


## ----fit_latent_nngp, fig.align="center", message=FALSE-----------------------
library(spNNGP)

n.samples <- 500

starting <- list("phi"=phi, "sigma.sq"=5, "tau.sq"=1)

tuning <- list("phi"=0.5, "sigma.sq"=0.5, "tau.sq"=0.5)

priors <- list("phi.Unif"=c(3/1, 3/0.01), "sigma.sq.IG"=c(2, 5), "tau.sq.IG"=c(2, 1))

cov.model <- "exponential"

m.s <- spNNGP(y~x-1, coords=coords, starting=starting, method="latent", n.neighbors=5,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, return.neighbor.info=TRUE, n.omp.threads=2)

round(summary(m.s$p.beta.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.s$p.theta.samples)$quantiles[,c(3,1,5)],2)
plot(w, apply(m.s$p.w.samples, 1, median),  xlab="True w", ylab="Posterior median w")


## ----look_at_neighbor_info----------------------------------------------------
names(m.s$neighbor.info)
ord <- m.s$neighbor.info$ord
n.indx <- m.s$neighbor.info$n.indx

s <- 50
plot(m.s$coords[ord,], xlab="Easting", ylab="Northing")
points(m.s$coords[ord,][s,,drop=FALSE], cex=2, col="red")
points(m.s$coords[ord,][n.indx[[s]],,drop=FALSE], pch=19, col="blue")
abline(v=m.s$coords[ord,][s,1], lty=3)

#for(s in 2:n){
#  plot(m.s$coords[ord,], xlab="Easting", ylab="Northing")
#  points(m.s$coords[ord,][s,,drop=FALSE], cex=2, col="red")
#  points(m.s$coords[ord,][n.indx[[s]],,drop=FALSE], pch=19, col="blue")
#  abline(v=m.s$coords[ord,][s,1], lty=3)
#  readline(prompt = "Pause. Press <Enter> to continue...")
#}



## ----harvard_forest_data, message=FALSE---------------------------------------
library(geoR)
library(raster)
library(leaflet)

load("CHM.rda")

CHM <- CHM[CHM[,3]>0,]

set.seed(1)
mod <- sample(1:nrow(CHM), 25000)
ho <- sample((1:nrow(CHM))[-mod], 10000)

CHM.mod <- CHM[mod,]
CHM.ho <- CHM[ho,]


## ----reproject_and_map, warning=FALSE-----------------------------------------
chm.r <- rasterFromXYZ(CHM)
proj4string(chm.r) <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
chm.r.ll <- projectRaster(chm.r, crs="+proj=longlat +datum=WGS84")

pal <- colorNumeric(rev(terrain.colors(50)), domain = values(chm.r.ll), na.color = "transparent")

base.map <- leaflet(width="100%") %>%
    addProviderTiles("Esri.WorldImagery", group="Satellite") %>%
    addProviderTiles("Esri.WorldShadedRelief", group="Terrain")

base.map %>%
    addRasterImage(chm.r.ll, colors = pal, opacity = 1, group="Canopy height") %>%
    addLegend("bottomright", pal = pal, values = values(chm.r.ll), opacity = 1, title = "<center>Canopy height (m)</center>") %>%
    addLayersControl(
        baseGroup = c("Satellite", "Terrain"),
        overlayGroups = c("Canopy height"),
        options = layersControlOptions(collapsed = FALSE)
    )


## ----variogram_eda, fig.align="center"----------------------------------------
sub <- 1:10000

#note, max intersite distance is ~1.5km
v <- variog(coords=CHM.mod[sub,1:2], data=CHM.mod[sub,3], uvec=(seq(0, 500, length=30))) 

plot(v, xlab="Distance (m)")


## ----fit_response_nngp--------------------------------------------------------
n.samples <- 1000

starting <- list("phi"=3/50, "sigma.sq"=15, "tau.sq"=2.5)

tuning <- list("phi"=0.05, "sigma.sq"=0.01, "tau.sq"=0.01)

priors <- list("phi.Unif"=c(3/1000, 3/10), "sigma.sq.IG"=c(2, 10), "tau.sq.IG"=c(2, 5))

cov.model <- "exponential"

##Response model 
m.r <- spNNGP(CHM.mod[,3] ~ 1, coords=CHM.mod[,1:2], starting=starting, method="response", n.neighbors=10,
              tuning=tuning, priors=priors, cov.model=cov.model,
              n.samples=n.samples, n.omp.threads=2, n.report=100)

round(summary(m.r$p.beta.samples)$quantiles[c(3,1,5)],2)
round(summary(m.r$p.theta.samples)$quantiles[,c(3,1,5)],3)

plot(m.r$p.beta.samples)
plot(m.r$p.theta.samples)

m.r$run.time



## ----predict_for_holdout, fig.align="center"----------------------------------
burn.in <- floor(0.5*n.samples)

p.r <- predict(m.r, X.0 = as.matrix(rep(1,nrow(CHM.ho))), coords.0 = CHM.ho[,c("x","y")],
                 sub.sample=list(start=burn.in, thin=2), n.report=5000, n.omp.threads=2)

y.hat.r <- apply(p.r$p.y.0, 1, mean)

plot(CHM.ho[,3], y.hat.r, main="Response NNGP model", xlab="True canopy height", ylab="Posterior predictive distribution mean")

