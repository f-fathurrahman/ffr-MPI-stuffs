Pkg.activate("CLANG", shared=true)
using Clang.Generators

include_dir = normpath(".")

args = get_default_args()
push!(args, "-I$include_dir")

# only wrap libclang headers in include/clang-c
headers = ["isddft.h"]

# create context
ctx = create_context(headers, args)

# run generator
build!(ctx)