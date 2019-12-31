cc_library(
    name = "jni_headers",
    srcs = [
        "@bazel_tools//tools/jdk:jni_header",
        "@bazel_tools//tools/jdk:jni_md_header-linux",
    ],
    includes = [
        "external/bazel_tools/tools/jdk/include",
        "external/bazel_tools/tools/jdk/include/linux",
    ],
    linkstatic = 1,
    visibility = [
        "//visibility:public",
    ],
)
