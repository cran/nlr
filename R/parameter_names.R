parameter.names<-
function(formula, data)
{
  pn <- names(attr(data, "parameters"))
  if(!length(pn))
    #pn <- database.attr("parameters")
  vars <- all.vars(formula)
  dnames <- names(data)
  ok <- !match(vars, dnames, 0.)
  if(length(pn))
    ok <- ok | match(vars, pn, 0.) > 0.
  vars[ok]
}
