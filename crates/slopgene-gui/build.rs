//! Build script of the desktop app. On Windows it puts the app icon into slopgene-gui.exe, so
//! Explorer, the Start menu, shortcuts and the taskbar show it. The icon is cosmetic: without a
//! resource compiler, or when compiling the resource fails, the build goes on with a warning.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=slopgene-gui.rc");
    println!("cargo:rerun-if-changed=assets/icon/slopgene.ico");
    embed_icon();
}

/// Compile `slopgene-gui.rc` and link it into the program, when building for Windows on Windows.
/// The resource compiler is `rc.exe` from the Windows SDK with MSVC, and `windres` with MinGW.
#[cfg(windows)]
fn embed_icon() {
    // A build script runs on the machine doing the build, which may build for another system.
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }
    if let Err(e) = embed_resource::compile("slopgene-gui.rc", embed_resource::NONE).manifest_optional() {
        println!("cargo:warning=slopgene-gui.exe was built without its icon: {e}");
    }
}

/// Other machines have no Windows resource compiler (and no `embed-resource`), so they skip it.
#[cfg(not(windows))]
fn embed_icon() {}
