package(default_visibility = ["//test:__subpackages__"])

cc_binary(
    name = "ASTRID",
    srcs = [
        "astrid.cpp",
    ],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods:DistanceMethods",
        "@phylokit",
    ],
)

cc_library(
    name = "Args",
    srcs = ["help.hpp"],
    hdrs = ["Args.hpp"],
)

cc_library(
    name = "Multind",
    srcs = ["multind.cpp"],
    hdrs = ["multind.hpp"],
    deps = [
        "@boost//:tokenizer",
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_library(
    name = "Octal",
    srcs = ["octal.cpp"],
    hdrs = ["octal.hpp"],
    deps = [
        "@phylokit",
    ],
)
