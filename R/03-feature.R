## 03-features.R
## Copyright (C) 2022-4 Kevin R. Coombes, RB McGee, and Jake Reed
## LICENSE: Perl Artistic License 2.0

Feature <- function(values, name, colors, meaning, ...) {
  colRamp <- colorRamp2(range(values, na.rm = TRUE),
                        colors = colors, ...)
  new("Feature",
      name = name,
      values = values,
      colRamp = colRamp,
      meaning = meaning)
}

setMethod("plot", c("Feature", "matrix"), function(x, y, pch = 16, ...) {
  plot(y, ...)
  points(x, y, pch = pch)
})

setMethod("points", "Feature", function(x, view, ...) {
  points(view, col = x@colRamp(x@values), ...)
})

