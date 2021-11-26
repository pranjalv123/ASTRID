cc_library(
    name = "DistanceMethods",
    hdrs = ["DistanceMethods.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":FastMEInterface",
        ":UPGMA",
    ],
)

cc_library(
    name = "FastMEInterface",
    srcs = ["FastMEInterface.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        "//third_party/fastme",
        "@phylokit",
    ],
)

cc_library(
    name = "UPGMA",
    srcs = ["UPGMA.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        "@com_github_google_glog//:glog",
        "@phylokit",
    ],
)

cc_library(
    name = "whereami",
    srcs = ["whereami.c"],
    hdrs = ["whereami.h"],
    copts=["-std=c99"]
)
