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

config_setting(
    name = "is_windows",
    constraint_values = [
        "@platforms//os:windows",
    ],
)

config_setting(
    name = "is_linux",
    constraint_values = [
        "@platforms//os:linux",
    ],
)

config_setting(
    name = "is_mac",
    constraint_values = [
        "@platforms//os:osx",
    ],
)
