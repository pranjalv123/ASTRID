package(default_visibility = ["//test:__subpackages__"])

cc_binary(
    name = "ASTRID-phydstar",
    srcs = [
        "astrid.cpp",
        "@bazel_tools//tools/jdk",
    ],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods",
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_binary(
    name = "ASTRID",
    srcs = [
        "astrid.cpp",
    ],
    deps = [
        ":Args",
        ":Multind",
        ":Octal",
        "//src/DistanceMethods:nophydstar",
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
