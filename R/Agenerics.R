register_generic = \(name, def = NULL) {
  if (!methods::isGeneric(name)){
    if (is.null(def)) {
      def = eval(bquote(function(data, ...) standardGeneric(.(name))))
    }
    methods::setGeneric(name, def)
  }
}

for (gen in c("embedded", "simplex", "smap",
              "ccm", "pcm", "multispatialccm")) {
  register_generic(gen)
}
