cc_library(
    name = "DistanceMethods",
    hdrs = ["DistanceMethods.hpp"],
    copts = ["-DUSE_PHYDSTAR"],
    visibility = ["//visibility:public"],
    deps = [
        "PhyDStar",
        ":FastMEInterface",
        ":RapidNJInterface",
        ":UPGMA",
    ],
)

cc_library(
    name = "nophydstar",
    hdrs = ["DistanceMethods.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "PhyDStar",
        ":FastMEInterface",
        ":RapidNJInterface",
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
    name = "RapidNJInterface",
    srcs = ["RapidNJInterface.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        "//third_party/rapidNJ",
        "@phylokit",
    ],
)

cc_library(
    name = "PhyDStar",
    srcs = [
        "phydstar.cpp",
    ],
    hdrs = ["DistanceMethods.hpp"],
    deps = [
        ":whereami",
        "//:jni_headers",
        "@phylokit",
    ],
)

cc_library(
    name = "UPGMA",
    srcs = ["UPGMA.cpp"],
    hdrs = ["DistanceMethods.hpp"],
    deps = ["@phylokit"],
)

cc_library(
    name = "whereami",
    srcs = ["whereami.c"],
    hdrs = ["whereami.h"],
)