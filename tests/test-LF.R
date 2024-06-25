library(RPointCloud)
data(CLL)
view <- cmdscale(daisydist)
circular <- angleMeans(view, ripdiag, NULL, clinical)
lf <- LoopFeature(circular)
round(sort(lf@Kstat), 4)
opar <- par(mai = c(0.82, 0.2, 0.82, 1.82))
image(lf, main = "Clinical Factors in CLL")
par(opar)

plot(lf, "mutation.status")
opar <- par(mfcol = c(2,2))
plot(lf, c("Serum.beta.2.microglobulin", "Hemoglobin",
           "Rai.Stage", "LogLDH"))
par(opar)

