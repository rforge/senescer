\name{senesceR-package}
\alias{senesceR-package}
\alias{senesceR}
\docType{package}
\title{
Simulation of an in vitro population of cells undergoing stochastic transition from replication to senescence.
}
\description{
Simulates inoculating a small population of cells into a flask of limited capacity.  Cells divide stochastically until flask becomes confluent, whereupon a sample of cells is transferred to a fresh flask and division can continue.  Cells can divide according to the commitment theory of cellular senescence (Kirkwood and Holliday) or just have fixed (or distributed) division potentials.
}
\details{
\tabular{ll}{
Package: \tab senesceR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-03\cr
License: \tab Artistic-2.0\cr
}
The main function is simsen, which simulates the replication and senescence of a population of cells.  This function writes a tab-delimited text file.  Auxiliary functions for reading in results, analysing and plotting are also included (e.g. histPD, readfiles, makeplot).
}
\author{
Conor Lawless

Maintainer: Conor Lawless <conor.lawless@ncl.ac.uk>
}
\references{

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Carry out some simulations
for (i in 1:100) simsen(Capacity=1000,Prat=0.25,Divrate=0.5,Deltarep=1.0,Pc=0.275,Meandiv=43)

# Read in the results
bigdat=readfiles()

# Generate plots and write to file
pdf("PopulationResults.pdf")

for (f in unique(bigdat$filename)){
	dat=bigdat[bigdat$filename==f,]
	makeplot(dat)
}
histPD(bigdat)

dev.off()
}
