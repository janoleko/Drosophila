pv = function (coef, time, degree, L = 24){
  if (degree == 1) {
    return(coef[, 1] + coef[, 2] * sin(2 * pi * time/L) + 
             coef[, 3] * cos(2 * pi * time/L))
  }
  else if (degree == 2) {
    return(coef[, 1] + coef[, 2] * sin(2 * pi * time/L) + 
             coef[, 3] * sin(2 * pi * time * 2/L) + 
             coef[, 4] * cos(2 * pi * time/L) + 
             coef[, 5] * cos(2 * pi * time * 2/L))
  }
  else if (degree == 3) {
    return(coef[, 1] + coef[, 2] * sin(2 * pi * time/L) + 
             coef[, 3] * sin(2 * pi * time * 2/L) + 
             coef[, 4] * sin(2 * pi * time * 3/L) + 
             coef[, 5] * cos(2 * pi * time/L) + 
             coef[, 6] * cos(2 * pi * time * 2/L) + 
             coef[, 7] * cos(2 * pi * time * 3/L))
  }
}