## 00-classes.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

setClass("Feature",
         slots = c(name = "character",
                   values = "numeric",
                   colRamp = "function",
                   meaning = "character" # or data.frame?
                   ))

setClass("Cycle",
         slots = c(index = "matrix",
                   dimension = "numeric",
                   color = "character"))

setClass("LoopCircos",
         contains = "Cycle",
         slots = c(angles = "matrix",
                   colors = "list")
         )

setClass("Projection",
         slots = c(phi = "numeric",
                   psi = "numeric",
                   displayphi = "numeric",
                   displaypsi = "numeric",
                   value = "matrix",
                   feature = "Feature"))

