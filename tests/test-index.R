library(RPointCloud)
data(cytof)

fixed <- disentangle(Arip, AML10.node287)
identical(fixed, AML10.node287.rips)
all.equal(fixed, AML10.node287.rips)
identical(fixed$diagram, AML10.node287.rips$diagram)
X <-fixed$diagram[]
Y <- AML10.node287.rips$diagram[]
all(X==Y)
identical(X, Y)
all.equal(X, Y)
identical(fixed$birthLocation, AML10.node287.rips$birthLocation)
identical(fixed$deathLocation, AML10.node287.rips$deathLocation)
identical(fixed$cycleLocation, AML10.node287.rips$cycleLocation)
all.equal(fixed$cycleLocation, AML10.node287.rips$cycleLocation)

