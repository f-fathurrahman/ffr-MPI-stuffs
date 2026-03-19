Pkg.activate("CLANG", shared=true)
using Clang.Generators

include_dir = normpath("../src/include")

args = get_default_args()
push!(args, "-I$include_dir")

headers = detect_headers(include_dir, args)

options = load_options("generator.toml")
ctx = create_context(headers, args, options)

# will will print out Julia code
build!(ctx)