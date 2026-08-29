fn main() {
    let zstd = pkg_config::Config::new()
        .probe("libzstd")
        .expect("libzstd is required for V10");
    let curl = pkg_config::Config::new()
        .probe("libcurl")
        .expect("libcurl is required for V10");
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++14")
        .file("src/v10_bridge.cpp")
        .file("../straw/C++/straw_v10.cpp")
        .include("../straw/C++");
    for path in zstd.include_paths.iter().chain(curl.include_paths.iter()) {
        build.include(path);
    }
    build.compile("straw_v10_bridge");
    println!("cargo:rerun-if-changed=src/v10_bridge.cpp");
    println!("cargo:rerun-if-changed=../straw/C++/straw_v10.cpp");
    println!("cargo:rerun-if-changed=../straw/C++/straw_v10.h");
    println!("cargo:rerun-if-changed=../straw/C++/v10_binary.h");
}
