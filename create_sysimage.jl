using PackageCompiler
create_sysimage(sysimage_path="./binder/sysimage_path",
                precompile_execution_file = "warmup.jl",
                replace_default = true,
                cpu_target = PackageCompiler.default_app_cpu_target())

