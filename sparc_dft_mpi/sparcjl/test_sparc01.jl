const libsparc = "./libsparc.so"

include("sparc_structs.jl")

function Initialize(pSPARC, argc, argv)
    ccall((:Initialize, libsparc), Cvoid, (Ptr{SPARC_OBJ}, Cint, Ptr{Ptr{Cchar}}), pSPARC, argc, argv)
end

