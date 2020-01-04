load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "9f9fb8b2f0213989247c9d5c0e814a8451d18d7f",
    remote = "https://github.com/nelhage/rules_boost",
    shallow_since = "1570056263 -0700",
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()

git_repository(
    name = "com_github_gflags_gflags",
    commit = "e171aa2d15ed9eb17054558e0b3a6a413bb01067",
    remote = "https://github.com/gflags/gflags.git",
    shallow_since = "1541971260 +0000",
)

git_repository(
    name = "com_github_google_glog",
    commit = "445af7ef7e6224340940c9b7d8e27c3b90d00101",
    remote = "https://github.com/google/glog/",
    shallow_since = "1553223106 +0900",
)

git_repository(
    name = "catch2",
    commit = "c5538476052dfe9d3ff2325198b1a8f42fc10669",
    remote = "https://github.com/evanmoran/catch2-bazel/",
    shallow_since = "1530732979 -0700",
)

git_repository(
    name = "phylokit",
    commit = "aac94d561d24cb747a1e8381ffa2da0a0a0d72b6",
    remote = "https://github.com/pranjalv123/phylokit/",
)
