# io
#  * give prefixes to the output files

out <- function(f){
  if(!dir.exists("output")){
    dir.create("output")
  }
  file.path("output", f)
}
