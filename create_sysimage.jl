using PackageCompiler
create_sysimage(["Pluto", "PlutoUI", "Graphs", "GraphPlot", "Compose", "QuantumAlgebra", "MetaGraphsNext"];
                sysimage_path="sysimage.so",
                precompile_execution_file = "warmup.jl",
                cpu_target = PackageCompiler.default_app_cpu_target())

