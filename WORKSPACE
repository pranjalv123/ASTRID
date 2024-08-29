load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "8fa193c4e21daaa2d46ff6b9c2b5a2de70b6caa1",
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
    commit = "1aba1a93495077267e1db1d138938d9f80678c98",
    remote = "https://github.com/pranjalv123/phylokit/",
)

git_repository(
    name = "platforms",
    commit = "d5fcb55bf8201539d0c3cb3a882c475330983156",
    remote = "https://github.com/bazelbuild/platforms.git",
)
